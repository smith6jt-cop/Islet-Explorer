// Native interpretation of regenerated cores and their exact centroid memberships.
// Args: candidate GeoJSON, native centroids TSV, case ID, fresh output directory.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.roi.PolygonROI
import qupath.lib.roi.GeometryROI
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.MultiPolygon
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.index.strtree.STRtree
import java.security.MessageDigest
import java.util.zip.GZIPOutputStream

def argv=this.binding.getVariable('args') as List
if(argv.size()!=4) throw new IllegalArgumentException('Provide candidates, centroids, case ID, fresh output')
def candidate=new File(argv[0].toString()).canonicalFile
def centroidFile=new File(argv[1].toString()).canonicalFile
def cid=argv[2].toString();def out=new File(argv[3].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve previous exact geometry export')
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def cores=PathIO.readObjects(candidate).sort{it.getName()}
if(cores.isEmpty() || cores.any{!it.getROI().getGeometry().isValid()}) throw new IllegalStateException('Missing or invalid regenerated cores')
if(cores.collect{it.getName()}.toSet().size()!=cores.size()) throw new IllegalStateException('Duplicate core names')
def polyCoords={g ->
    def rings=[g.getExteriorRing().getCoordinates().collect{[it.x,it.y]}]
    for(int i=0;i<g.getNumInteriorRing();i++) rings.add(g.getInteriorRingN(i).getCoordinates().collect{[it.x,it.y]})
    rings
}
def gson=GsonTools.getInstance(true)
def features=cores.collect{object ->
    def feature=GsonTools.getInstance().fromJson(GsonTools.getInstance().toJson(object),Map.class)
    def roi=object.getROI()
    if(roi instanceof PolygonROI) {
        def coords=roi.getAllPoints().collect{[it.getX(),it.getY()]};coords.add(coords[0])
        feature.geometry=[type:'Polygon',coordinates:[coords]]
    } else if(roi instanceof GeometryROI) {
        def g=roi.getGeometry()
        if(g instanceof Polygon) feature.geometry=[type:'Polygon',coordinates:polyCoords(g)]
        else if(g instanceof MultiPolygon) feature.geometry=[type:'MultiPolygon',coordinates:(0..<g.getNumGeometries()).collect{polyCoords(g.getGeometryN(it))}]
        else throw new IllegalStateException('Unexpected area geometry')
    } else throw new IllegalStateException('Unsupported native ROI '+roi.getClass().getName())
    feature.get('properties').put('classification',[name:'Islet',color:[0,255,255]])
    feature.get('properties').put('native_roi_type',roi.getClass().getSimpleName())
    feature
}
out.mkdirs()
new File(out,'cores.geojson.gz').withOutputStream{s ->new GZIPOutputStream(s).withWriter('UTF-8'){w ->
    GsonTools.getInstance().toJson([type:'FeatureCollection',features:features],w)
}}
def index=new STRtree();cores.each{index.insert(it.getROI().getGeometry().envelopeInternal,it)};index.build()
int count=0,coreCount=0
new File(out,'native_core_memberships.csv').withWriter('UTF-8'){w ->
    w.write('case_id,object_id,islet_key\n')
    centroidFile.withReader('UTF-8'){r ->
        def header=r.readLine().split('\t');def line
        while((line=r.readLine())!=null){
            def fields=line.split('\t');double x=Double.parseDouble(fields[1]),y=Double.parseDouble(fields[2]);count++
            def matches=index.query(new Envelope(x,x,y,y)).findAll{it.getROI().contains(x,y)}
            if(matches.size()>1) throw new IllegalStateException('Multiple core containment for '+fields[0])
            if(!matches.isEmpty()){coreCount++;w.write(cid+','+cid+':'+fields[0]+','+matches[0].getName()+'\n')}
        }
    }
}
new File(out,'manifest.json').text=gson.toJson([
    status:'native_geometry_and_memberships_exported',case_id:cid,source_candidates_sha256:sha(candidate),
    source_centroids_sha256:sha(centroidFile),cores:cores.size(),cells:count,core_cells:coreCount,
    cores_sha256:sha(new File(out,'cores.geojson.gz')),memberships_sha256:sha(new File(out,'native_core_memberships.csv')),
    source_project_modified:false,qupath_version:qupath.lib.common.GeneralTools.getVersion(),
    roi_types:cores.countBy{it.getROI().getClass().getSimpleName()},
    geometry_precision:'Native PolygonROI vertices or exact GeometryROI coordinates; no additional JTS rounding'])
println "${cid}: ${cores.size()} regenerated native cores; ${coreCount} core cells"
