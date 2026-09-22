// Exercise the same morphology on the original 15 reference geometries.
// This is a geometry/transfer check, not new classifier inference or ground truth.
// Args: exported reference root, fresh output root, geometry helper.
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import java.security.MessageDigest

def argv=this.binding.getVariable('args') as List
if(argv.size()!=3) throw new IllegalArgumentException('Provide reference root, fresh output root, geometry helper')
def input=new File(argv[0].toString()).canonicalFile
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve previous reference checks')
def helper=new File(argv[2].toString()).canonicalFile
def rules=new GroovyShell(getClass().getClassLoader()).evaluate(helper)
def gson=GsonTools.getInstance(true)
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def refs=gson.fromJson(new File(input,'reference_manifest.json').text,Map.class)
def expected=['112','6356','6450','6479','6505','6516','6521','6533','6534','6538','6548','6549','6550','6551','6563'].toSet()
if(refs.donors.size()!=15 || refs.donors.collect{it.case_id}.toSet()!=expected)
    throw new IllegalStateException('Expected exactly the original 15 donors')
out.mkdirs()
def summaries=[]
refs.donors.sort{a,b ->a.case_id.toInteger()<=>b.case_id.toInteger()}.each{ref ->
    def source=new File(ref.reference_annotations.toString())
    if(sha(source)!=ref.reference_sha256) throw new IllegalStateException('Reference annotations changed')
    def objects=PathIO.readObjects(source)
    def islets=objects.findAll{it.getPathClass()?.toString()=='Islet'}
    def tissue=GeometryTools.union(objects.findAll{it.getPathClass()?.toString()=='Tissue'}.collect{it.getROI().getGeometry()})
    def result=rules.refine(islets.collect{it.getROI().getGeometry()},tissue,
        ref.image_shape_cyx[2].intValue(),ref.image_shape_cyx[1].intValue(),ref.pixel_size_um.doubleValue(),ref.pixel_size_um.doubleValue())
    def folder=new File(out,ref.case_id.toString());folder.mkdirs()
    def revised=result.records.withIndex().collect{r,i ->
        def obj=PathObjects.createAnnotationObject(GeometryTools.geometryToROI(r.geometry,ImagePlane.getDefaultPlane()),PathClass.fromString('Adaptive reference geometry'))
        obj.setName(String.format('Reference_%s_%05d',ref.case_id,i+1))
        obj.getMetadata().put('Size category',r.size_category.toString())
        obj.getMeasurementList().put('Parent index',r.parent_index as double)
        obj
    }
    def target=new File(folder,'reference_adaptive.geojson')
    PathIO.exportObjectsAsGeoJSON(target,revised,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
    new File(folder,'parent_audit.json').text=gson.toJson(result.audit)
    def summary=[case_id:ref.case_id,source_reference:source.canonicalPath,
        source_reference_sha256:sha(source),source_reference_unchanged:true,
        output:target.canonicalPath,output_sha256:sha(target),summary:result.summary]
    new File(folder,'verification.json').text=gson.toJson(summary)
    summaries.add(summary)
    println 'Checked original donor '+ref.case_id+': '+gson.toJson(result.summary)
}
new File(out,'verification.json').text=gson.toJson([status:'passed',donor_count:summaries.size(),
    geometry_helper_sha256:sha(helper),config:rules.config,donors:summaries,
    scope:'Morphology on preserved reference shapes; no image classifier refitting, inference or biological accuracy estimate',
    source_native_files_opened:false,source_project_modified:false,accepted_islets:0])
