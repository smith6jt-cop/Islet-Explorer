// Promote only the explicitly approved, hash-bound HPAP review in an independent copy.
// Args: staged six-donor project root, fresh HPAP export folder.
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import qupath.lib.analysis.DistanceTools
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.union.UnaryUnionOp
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.nio.charset.StandardCharsets
import java.security.MessageDigest
import java.util.zip.GZIPOutputStream

def argv=this.binding.getVariable('args') as List
if(argv.size()!=2) throw new IllegalArgumentException('Provide staged project and fresh export folder')
def root=new File(argv[0].toString()).canonicalFile
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve existing promotion/export evidence')
def gson=GsonTools.getInstance(true)
def approval=gson.fromJson(new File(root,'boundary_approval.json').text,Map.class)
def staging=gson.fromJson(new File(root,'staging_manifest.json').text,Map.class)
if(approval.case_id!='78' || approval.decision!='owner_visually_approved' || !approval.user_statement)
    throw new IllegalStateException('Explicit HPAP boundary approval required')
def donor=staging.donors.find{it.case_id=='78'}
def source=new File(approval.source_review_data.toString()).canonicalFile
def target=new File(donor.review_data.toString()).canonicalFile
def sha={File file ->
    def d=MessageDigest.getInstance('SHA-256')
    file.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
if(!target.toPath().startsWith(root.toPath()) || Files.isSameFile(source.toPath(),target.toPath()))
    throw new IllegalStateException('An independent writable project copy is required')
if(sha(source)!=approval.source_review_data_sha256 || sha(target)!=approval.source_review_data_sha256)
    throw new IllegalStateException('Reviewed or staged native data changed')
def changingDistances=['Signed distance to annotation Islet µm','Signed distance to annotation IsletExpanded µm'].toSet()
def hashCells={hierarchy ->
    def d=MessageDigest.getInstance('SHA-256')
    hierarchy.getCellObjects().sort{it.getID().toString()}.each{cell ->
        def measurements=new TreeMap(cell.getMeasurements())
        changingDistances.each{measurements.remove(it)}
        def r=[id:cell.getID().toString(),classification:cell.getPathClass()?.toString(),
            roi:cell.getROI().getGeometry().toText(),nucleus:cell.getNucleusROI()?.getGeometry()?.toText(),measurements:measurements]
        d.update(GsonTools.getInstance().toJson(r).getBytes(StandardCharsets.UTF_8));d.update((byte)10)
    }
    d.digest().encodeHex().toString()
}
println 'Reading the exact owner-reviewed HPAP native hierarchy'
def imageData=PathIO.readImageData(target)
def hierarchy=imageData.getHierarchy()
def cells=hierarchy.getCellObjects().toList()
def originals=hierarchy.getAnnotationObjects().collectEntries{[(it.getID().toString()):it.getROI().getGeometry().toText()]}
def cores=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString()==approval.source_review_object_class}.sort{it.getName()}
if(cores.isEmpty() || cores.any{!it.getROI().isArea() || !it.getROI().getGeometry().isValid()})
    throw new IllegalStateException('Reviewed cores are missing or invalid')
if(cores.collect{it.getName()}.toSet().size()!=cores.size()) throw new IllegalStateException('Duplicate reviewed core names')
if(cells.size()!=donor.cell_count.intValue()) throw new IllegalStateException('Raw cell count changed')
if(hierarchy.getAnnotationObjects().any{it.getPathClass()?.toString() in ['Islet','IsletExpanded']})
    throw new IllegalStateException('Previous active islet classes must already be archived in the reviewed project')
double sx=imageData.getServer().getPixelCalibration().getPixelWidthMicrons()
double sy=imageData.getServer().getPixelCalibration().getPixelHeightMicrons()
if(Math.abs(sx-.3774)>1e-9 || Math.abs(sy-.3774)>1e-9) throw new IllegalStateException('Unexpected spatial calibration')
def tissues=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString()=='Tissue'}
def tissue=UnaryUnionOp.union(tissues.collect{it.getROI().getGeometry()})
def excludes=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString()=='Exclude'}
if(!excludes.isEmpty()) tissue=tissue.difference(UnaryUnionOp.union(excludes.collect{it.getROI().getGeometry()}))
String before=hashCells(hierarchy)
out.mkdirs()
def originalParents=cells.collectEntries{[(it.getID().toString()):it.getParent()?.getID()?.toString()]}
def index=new STRtree()
cores.each{core ->index.insert(core.getROI().getGeometry().envelopeInternal,core)}
index.build()
cores.each{core ->
    core.setPathClass(PathClass.fromString('Islet'))
    core.setColor(0,255,255)
    core.getMetadata().put('Boundary review','Owner visually approved on 2026-09-18')
    core.getMetadata().put('Review status','Accepted boundary; phenotype gates remain independently unvalidated')
    core.getMetadata().put('Approval record SHA256',sha(new File(root,'boundary_approval.json')))
}
def groups=[:].withDefault{[]}
int coreCells=0,changedParents=0
new File(out,'native_core_memberships.csv').withWriter('UTF-8'){writer ->
    writer.writeLine('case_id,object_uuid,cell_id,islet_key,islet_uuid')
    cells.each{cell ->
        double x=cell.getROI().getCentroidX(),y=cell.getROI().getCentroidY()
        def matches=index.query(new Envelope(x,x,y,y)).findAll{it.getROI().contains(x,y)}
        if(matches.size()>1) throw new IllegalStateException('Cell belongs to more than one reviewed core: '+cell.getID())
        def parent=matches.isEmpty() ? tissues.find{it.getROI().contains(x,y)} : matches[0]
        if(parent==null) parent=hierarchy.getRootObject()
        if(parent.getID()?.toString()!=originalParents[cell.getID().toString()]) changedParents++
        parent.addChildObject(cell)
        if(!matches.isEmpty()) {
            coreCells++;groups[matches[0].getID().toString()].add(cell)
            writer.writeLine('78,'+cell.getID()+',78:'+cell.getID()+','+matches[0].getName()+','+matches[0].getID())
        }
    }
}
def markers=['CD31','CD8a','CD3e','SMA','CD4','CD34','HLADR','PDPN','panCK','ECAD','CD163','SST','ColIV','VIM','CD20','LGALS3',
    'B3TUBB','GCG','KRT14','GAP43','CD35','CHGA','PGP9.5','INS','CD44','CD45','NaKATPase','BCatenin','CD68','BActin','CK19']
cores.each{core ->
    def members=groups[core.getID().toString()]
    def measured=members.findAll{cell ->markers.every{Double.isFinite(cell.getMeasurementList().getOrDefault('Cell: '+it+': Mean',Double.NaN))}}
    core.getMeasurementList().put('Rebuilt core: raw cells',members.size())
    core.getMeasurementList().put('Rebuilt core: measured cells',measured.size())
    core.getMeasurementList().put('Rebuilt core: wholly unmeasured cells',members.size()-measured.size())
    core.getMeasurementList().put('Rebuilt core: section area um2',core.getROI().getArea()*sx*sy)
    markers.each{marker ->
        def values=measured.collect{Math.log1p(it.getMeasurementList().get('Cell: '+marker+': Mean'))}.sort()
        int n=values.size()
        double median=n==0 ? Double.NaN : (n%2 ? values[n.intdiv(2)] : (values[n.intdiv(2)-1]+values[n.intdiv(2)])/2d)
        core.getMeasurementList().put('Rebuilt core: median log1p '+marker,median)
    }
}
// Native context polygons are full core+20 um expansions. Regional analysis
// separately subtracts the union of all cores to obtain the peri-only band.
def expansions=cores.collect{core ->
    def g=core.getROI().getGeometry().buffer(20d/sx,16).intersection(tissue)
    def obj=PathObjects.createAnnotationObject(GeometryTools.geometryToROI(g,ImagePlane.getDefaultPlane()),PathClass.fromString('IsletExpanded'))
    obj.setName(core.getName()+'_expanded20')
    obj.setColor(220,190,60)
    obj.getMetadata().put('Core UUID',core.getID().toString())
    obj.getMetadata().put('Definition','Tissue-clipped core + 20 um; subtract all cores for peri-only analysis')
    obj
}
hierarchy.getRootObject().addChildObjects(expansions)
println "Reparented ${coreCells} core cells; rebuilding signed distances for ${cells.size()} cells"
DistanceTools.centroidToBoundsDistance2D(cells,cores,sx,sy,'Signed distance to annotation Islet µm',true)
DistanceTools.centroidToBoundsDistance2D(cells,expansions,sx,sy,'Signed distance to annotation IsletExpanded µm',true)
hierarchy.fireHierarchyChangedEvent(this)
def activeAnnotations=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString() in ['Islet','IsletExpanded','Tissue','Exclude','Nerve','Capillary','Lymphatic']}
new File(out,'annotations.geojson.gz').withOutputStream{stream ->
    new GZIPOutputStream(stream).withWriter('UTF-8'){writer ->
        GsonTools.getInstance().toJson([type:'FeatureCollection',features:activeAnnotations],writer)
    }
}
PathIO.exportObjectsAsGeoJSON(new File(out,'accepted_islets.geojson'),cores,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
imageData.setProperty('Endocrine pixel review','HPAP boundaries owner-approved 2026-09-18; native core cell parents and signed distances rebuilt. Phenotype validation remains separate.')
def temporary=new File(target.parentFile,'approved_boundaries.qpdata')
PathIO.writeImageData(temporary,imageData)
imageData.close();imageData=null;hierarchy=null;cells=null;cores=null;expansions=null;activeAnnotations=null;groups=null;System.gc()
def saved=PathIO.readHierarchy(temporary)
String after=hashCells(saved)
if(before!=after || saved.getCellObjects().size()!=donor.cell_count.intValue())
    throw new IllegalStateException('Cell shapes, raw measurements, classifications or unrelated distances changed')
def annotations=saved.getAnnotationObjects().collectEntries{[(it.getID().toString()):it]}
originals.each{id,wkt ->if(annotations[id]?.getROI()?.getGeometry()?.toText()!=wkt)throw new IllegalStateException('An existing annotation shape changed')}
if(saved.getCellObjects().count{it.getParent()?.getPathClass()?.toString()=='Islet'}!=coreCells)
    throw new IllegalStateException('Core parent assignments changed on serialization')
if(sha(source)!=approval.source_review_data_sha256) throw new IllegalStateException('Reviewed source changed during promotion')
Files.move(temporary.toPath(),target.toPath(),StandardCopyOption.REPLACE_EXISTING)
def summaryFile=new File(target.parentFile,'summary.json')
def summary=gson.fromJson(summaryFile.text,Map.class)
summary.hierarchy.annotationClassificationCounts=saved.getAnnotationObjects().countBy{it.getPathClass()?.toString() ?: 'Unclassified'}
summary.hierarchy.objectTypeCounts.Annotations=saved.getAnnotationObjects().size()
summary.hierarchy.nObjects=saved.nObjects()+1;summary.timestamp=System.currentTimeMillis();summaryFile.text=gson.toJson(summary)
def projectFile=new File(root,'project.qpproj');def project=gson.fromJson(projectFile.text,Map.class)
def entry=project.images.find{it.entryID.intValue()==42}
entry.metadata['Boundary review']='Owner visually approved 2026-09-18; exact reviewed shapes promoted'
entry.metadata['Delivery status']='HPAP boundaries accepted; native cell membership rebuilt; phenotype validation remains separate'
entry.metadata['Measurements']='Core membership, 31 core marker medians and islet distances rebuilt; regional tables generated separately'
projectFile.text=gson.toJson(project)
def evidence=[status:'approved_boundaries_promoted_and_native_roundtrip_passed',case_id:'78',
    accepted_islets:summary.hierarchy.annotationClassificationCounts.Islet,cells:saved.getCellObjects().size(),
    core_cells:coreCells,changed_cell_parents:changedParents,
    unchanged_cell_content_sha256_before:before,unchanged_cell_content_sha256_after:after,
    raw_cell_measurements_and_geometry_preserved:true,all_original_annotation_geometry_preserved:true,
    recalculated_cell_measurements:changingDistances,cell_parents_rebuilt:true,
    source_review_unmodified:true,source_data_sha256:approval.source_review_data_sha256,
    output_native_sha256:sha(target),output_project_sha256:sha(projectFile),
    core_geometry_sha256:sha(new File(out,'accepted_islets.geojson')),
    geometry_sha256:sha(new File(out,'annotations.geojson.gz')),approval_sha256:sha(new File(root,'boundary_approval.json')),
    phenotype_gates_validated_by_this_approval:false]
new File(out,'promotion_manifest.json').text=gson.toJson(evidence)
println gson.toJson(evidence)
