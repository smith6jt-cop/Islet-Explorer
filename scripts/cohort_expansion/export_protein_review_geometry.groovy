// Read-only snapshot of the six-donor review. Does not approve candidates.
// Args: current review project directory, fresh output directory.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.roi.PolygonROI
import java.security.MessageDigest
import java.util.zip.GZIPOutputStream

def argv=this.binding.getVariable('args') as List
if(argv.size()!=2) throw new IllegalArgumentException('Provide review project directory and fresh output')
def project=new File(argv[0].toString()).canonicalFile
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve previous snapshots')
def gson=GsonTools.getInstance(true)
def delivery=gson.fromJson(new File(project,'delivery_manifest.json').text,Map.class)
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
out.mkdirs()
def records=[]
delivery.donors.each{donor ->
    def cid=donor.case_id.toString()
    def source=new File(donor.review_data.toString())
    def before=sha(source)
    def hierarchy=PathIO.readHierarchy(source)
    def wanted=cid=='78' ? 'Islet' : 'Endocrine pixel candidate'
    def cores=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString()==wanted}.sort{it.getID().toString()}
    if(cores.isEmpty()) throw new IllegalStateException('No review boundaries for '+cid)
    def aliases=[:]
    cores.eachWithIndex{obj,i ->aliases[obj.getID().toString()]='Islet_'+(i+1)}
    def annotations=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString() in ['Tissue','Exclude','Nerve','Capillary','Lymphatic']}
    def features=(cores+annotations).collect{obj ->
        def feature=GsonTools.getInstance().fromJson(GsonTools.getInstance().toJson(obj),Map.class)
        def roi=obj.getROI()
        if(roi instanceof PolygonROI){
            def coords=roi.getAllPoints().collect{[it.getX(),it.getY()]};coords.add(coords[0])
            feature.put('geometry',[type:'Polygon',coordinates:[coords]])
        } else {
            // JTS geometry serialization preserves doubles without QuPath's export rounding.
            def g=roi.getGeometry()
            def poly={p ->
                def rings=[p.getExteriorRing().getCoordinates().collect{[it.x,it.y]}]
                for(int i=0;i<p.getNumInteriorRing();i++)rings.add(p.getInteriorRingN(i).getCoordinates().collect{[it.x,it.y]})
                rings
            }
            if(g.getGeometryType()=='Polygon')feature.put('geometry',[type:'Polygon',coordinates:poly(g)])
            else if(g.getGeometryType()=='MultiPolygon')feature.put('geometry',[type:'MultiPolygon',coordinates:(0..<g.getNumGeometries()).collect{poly(g.getGeometryN(it))}])
            else throw new IllegalStateException('Non-polygon annotation '+g.getGeometryType())
        }
        if(obj.getPathClass()?.toString()==wanted) {
            def alias=aliases.get(obj.getID().toString())
            if(alias==null)throw new IllegalStateException('Missing stable boundary alias')
            def properties=feature.get('properties')
            properties.put('source_name',obj.getName())
            properties.put('source_object_id',obj.getID().toString())
            properties.put('name',alias)
            properties.put('classification',[name:'Islet',color:[0,255,255]])
            properties.put('native_roi_type',roi.getClass().getSimpleName())
            properties.put('review_status',cid=='78' ? 'owner_approved' : 'awaiting_owner_review')
        }
        feature
    }
    def dir=new File(out,cid);dir.mkdirs()
    def geometry=new File(dir,'annotations.geojson.gz')
    geometry.withOutputStream{s ->new GZIPOutputStream(s).withWriter('UTF-8'){w ->gson.toJson([type:'FeatureCollection',features:features],w)}}
    def centroid=new File(dir,'native_centroids.tsv')
    centroid.withWriter('UTF-8'){w ->
        w.write('object_id\tx_source_px\ty_source_px\n')
        hierarchy.getCellObjects().each{cell ->w.write([cell.getID().toString(),cell.getROI().getCentroidX(),cell.getROI().getCentroidY()].join('\t')+'\n')}
    }
    if(sha(source)!=before)throw new IllegalStateException('Native source changed during snapshot')
    records.add([case_id:cid,native_data:source.canonicalPath,native_sha256:before,
        geometry:geometry.canonicalPath,geometry_sha256:sha(geometry),
        centroids:centroid.canonicalPath,centroids_sha256:sha(centroid),
        cores:cores.size(),cells:hierarchy.getCellObjects().size(),
        boundary_status:cid=='78' ? 'owner_approved' : 'awaiting_owner_review',
        pixel_size_um:0.3774,source_modified:false,
        identity_mapping:'Islet_N assigned by sorted immutable native object UUID; source_name and source_object_id retained'])
    println "${cid}: ${cores.size()} boundaries snapshotted; source unchanged"
}
new File(out,'manifest.json').text=gson.toJson([status:'review_snapshot_not_approval',donors:records])
