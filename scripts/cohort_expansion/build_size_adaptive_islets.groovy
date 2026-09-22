// Apply IO60 morphology to saved endocrine foreground, using native QuPath/JTS.
// Args: prediction folder, FRESH output folder, geometry helper, IO60 reference script,
// optional reviewed technical exclusion FeatureCollection in source pixels.
// The forest and cutoff are fixed inputs. No source project is loaded or saved.
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.Geometry
import java.nio.file.Files
import java.security.MessageDigest

def argv = this.binding.getVariable('args') as List
if (!(argv.size() in [4,5])) throw new IllegalArgumentException('Provide prediction folder, fresh output folder, geometry helper, IO60 reference script, optional technical exclusions')
def input = new File(argv[0].toString()).canonicalFile
def out = new File(argv[1].toString()).canonicalFile
def helperFile = new File(argv[2].toString()).canonicalFile
def referenceFile = new File(argv[3].toString()).canonicalFile
if (out.exists()) throw new IllegalStateException('Preserve existing results and user review edits')
def gson = GsonTools.getInstance(true)
def sha = {File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def manifestFile = new File(input,'prediction_manifest.json')
def manifest = gson.fromJson(manifestFile.text,Map.class)
if (!(manifest.status in ['refined_pixel_prediction_unreviewed','unanimous_vote_sensitivity_not_selected_primary']))
    throw new IllegalStateException('Expected a documented endocrine-only prediction')
if (manifest.markers.toSet()!=['INS','GCG','SST'].toSet() || manifest.markers.size()!=3)
    throw new IllegalStateException('Unexpected foreground marker panel')
def vectors = gson.fromJson(new File(input,'vectorization_manifest.json').text,Map.class)
def rawFile = new File(input,'raw_pixel_regions.geojson')
def rawRecord = vectors.outputs.find{it.kind=='pixel_classes'}
if (vectors.source_prediction_manifest_sha256!=sha(manifestFile) || rawRecord==null || sha(rawFile)!=rawRecord.sha256)
    throw new IllegalStateException('Prediction and raw vector provenance differ')
def specFile = new File(manifest.source_specification.toString()).canonicalFile
if (sha(specFile)!=manifest.source_specification_sha256) throw new IllegalStateException('Source specification changed')
def spec = gson.fromJson(specFile.text,Map.class)
def tissueFile = new File(specFile.parentFile,'valid_tissue.geojson')
if (sha(tissueFile)!=spec.valid_tissue_sha256 || sha(tissueFile)!=vectors.tissue_sha256)
    throw new IllegalStateException('Tissue geometry changed')
def image = new File(spec.source_image.toString()).canonicalFile
if (image.canonicalPath!=new File(manifest.source_image.toString()).canonicalPath ||
    image.length()!=manifest.source_bytes.longValue() || image.lastModified()!=Math.floor(manifest.source_mtime_ns.doubleValue()/1e6).longValue())
    throw new IllegalStateException('Source image identity or fingerprint changed')
double sx=manifest.source_pixel_size_um.doubleValue(), sy=sx
if (sx!=spec.pixel_size_um.doubleValue() || spec.image_shape_cyx[2].intValue()!=manifest.source_width.intValue() ||
    spec.image_shape_cyx[1].intValue()!=manifest.source_height.intValue()) throw new IllegalStateException('Source spatial units differ')
// Explicitly check the referenced rule set; no dependence on the source
// project's hard-coded 0.5078 calibration or its CD99-dependent Islets2 model.
def expected = [edgeMarginUm:5d,SMALL_ISLET_AREA_UM2:2000d,MEDIUM_ISLET_AREA_UM2:10000d,
    SMALL_EROSION_UM:1d,MEDIUM_EROSION_UM:2d,LARGE_EROSION_UM:3d,DILATION_FACTOR:0.85d]
expected.each{name,value ->
    def match = referenceFile.text =~ /double\s+${name}\s*=\s*([\d.]+)/
    if (!match.find() || Double.parseDouble(match.group(1))!=value)
        throw new IllegalStateException('Reference morphology changed: '+name)
}
if (!referenceFile.text.contains('"minFragmentSizeMicrons":50.0,"maxHoleSizeMicrons":-1.0'))
    throw new IllegalStateException('Reference fragment/hole policy changed')
def rules = new GroovyShell(getClass().getClassLoader()).evaluate(helperFile)
def rawObjects = PathIO.readObjects(rawFile)
if (rawObjects.size()!=rawRecord.objects.intValue()) throw new IllegalStateException('Raw foreground count differs')
def tissue = gson.fromJson(tissueFile.text,Geometry.class)
def exclusionRecord=null
if(argv.size()==5) {
    def exclusionsFile=new File(argv[4].toString()).canonicalFile
    def exclusionObjects=PathIO.readObjects(exclusionsFile)
    if(exclusionObjects.any{it.getPathClass()?.toString()!='Exclude' || !it.getROI().getGeometry().isValid()})
        throw new IllegalStateException('Technical exclusions require valid Exclude annotations')
    def exclusionUnion=org.locationtech.jts.operation.union.UnaryUnionOp.union(exclusionObjects.collect{it.getROI().getGeometry()})
    double beforeArea=tissue.area
    if(exclusionUnion!=null && !exclusionUnion.isEmpty()) tissue=tissue.difference(exclusionUnion)
    if(tissue.isEmpty() || !tissue.isValid()) throw new IllegalStateException('Technical exclusion removed valid tissue')
    exclusionRecord=[path:exclusionsFile.canonicalPath,sha256:sha(exclusionsFile),
        excluded_tissue_area_um2:(beforeArea-tissue.area)*sx*sy,
        prediction_normalization_changed:false,role:'image artifact exclusion after fixed endocrine classification']
}
long started=System.nanoTime()
def result = rules.refine(rawObjects.collect{it.getROI().getGeometry()},tissue,
    manifest.source_width.intValue(),manifest.source_height.intValue(),sx,sy)
out.mkdirs()
new File(out,'assessable_tissue.geojson').text=gson.toJson(tissue)
def objects=[]
result.records.eachWithIndex{r,i ->
    def obj=PathObjects.createAnnotationObject(GeometryTools.geometryToROI(r.geometry,ImagePlane.getDefaultPlane()),
        PathClass.fromString('Endocrine pixel candidate'))
    obj.setName(String.format('AdaptiveIslet_%s_%05d',manifest.case_id,i+1))
    obj.setColor(0,255,255)
    ['area_um2','input_area_um2','erosion_um','dilation_um','foreground_area_um2',
        'filled_background_area_um2','foreground_fraction','parent_index'].each{key ->
        obj.getMeasurementList().put('Adaptive boundary: '+key,(r[key] as Number).doubleValue())
    }
    obj.getMetadata().put('Size category',r.size_category.toString())
    obj.getMetadata().put('Review status','Size-adaptive candidate; cell membership and regional measurements have not been recomputed')
    obj.getMetadata().put('Segmentation method','INS/GCG/SST foreground; IO60 size-adaptive erosion/dilation, 50 um2 fragments, enclosed-hole filling')
    obj.getMetadata().put('Model SHA256',manifest.model_sha256.toString())
    obj.getMetadata().put('Boundary method',rules.config.method.toString())
    objects.add(obj)
}
def overlay=new File(out,'review_candidates.geojson')
PathIO.exportObjectsAsGeoJSON(overlay,objects,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
def excluded=result.removed.collect{r ->
    def obj=PathObjects.createAnnotationObject(GeometryTools.geometryToROI(r.geometry,ImagePlane.getDefaultPlane()),PathClass.fromString('Endocrine fragment review'))
    obj.setName('Fragment_'+manifest.case_id+'_'+r.parent_index)
    obj.getMetadata().put('Reason',r.reason.toString());obj
}
PathIO.exportObjectsAsGeoJSON(new File(out,'excluded_fragments.geojson'),excluded,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
def auditFile=new File(out,'parent_audit.json');auditFile.text=gson.toJson(result.audit)
Files.copy(referenceFile.toPath(),new File(out,'reference_1_panc_objects.groovy').toPath())
Files.copy(helperFile.toPath(),new File(out,'size_adaptive_islet_geometry.groovy').toPath())
def parentStatus=manifest.status
manifest.status='size_adaptive_pixel_candidates_unreviewed'
manifest.parent_prediction_status=parentStatus
manifest.parent_prediction_manifest=manifestFile.canonicalPath
manifest.parent_prediction_manifest_sha256=sha(manifestFile)
manifest.raw_foreground_geojson=rawFile.canonicalPath
manifest.raw_foreground_sha256=sha(rawFile)
manifest.boundary_config=rules.config
manifest.technical_exclusions=exclusionRecord
manifest.assessable_tissue_path=new File(out,'assessable_tissue.geojson').canonicalPath
manifest.assessable_tissue_sha256=sha(new File(out,'assessable_tissue.geojson'))
manifest.boundary_policy=[active:rules.config.method,stage:'complete',raw_foreground_retained:true]
manifest.boundary_summary=result.summary
manifest.boundary_area_um2=result.summary.output_area_um2
manifest.boundary_connected_regions=objects.size()
manifest.boundary_method_reference=[path:referenceFile.canonicalPath,sha256:sha(referenceFile),
    helper_sha256:sha(helperFile),pixel_calibration_override_applied:false,
    reference_classifier_reused:false,reference_classifier_incompatibility:'Islets2 requires CD99, absent from this cohort',
    distance_conversion:'Source script buffers pixel coordinates directly; adaptation applies stated distances in micrometres']
manifest.parent_raster_sha256=manifest.remove('raster_sha256')
manifest.boundary_runtime_seconds=(System.nanoTime()-started)/1e9
manifest.boundary_review_required=true
manifest.operating_point=manifest.operating_point ?: [probability_cutoff:manifest.config.probability_cutoff,
    status:'unchanged parent cutoff; morphology does not validate foreground classification']
def newManifest=new File(out,'prediction_manifest.json');newManifest.text=gson.toJson(manifest)
new File(out,'vectorization_manifest.json').text=gson.toJson([
    status:'native_size_adaptive_refinement_complete',qupath_version:qupath.lib.common.GeneralTools.getVersion(),
    source_prediction_manifest_sha256:sha(newManifest),tissue_sha256:sha(tissueFile),
    prediction_downsample:manifest.prediction_downsample,
    outputs:[[kind:'boundary_regions',objects:objects.size(),geometry_valid:true,
        clipped_area_source_px2:result.summary.output_area_um2/(sx*sy),path:overlay.canonicalPath,sha256:sha(overlay)]],
    parent_audit_sha256:sha(auditFile),original_project_modified:false,source_native_objects_modified:false,accepted_islets:0])
println gson.toJson(result.summary)
println 'Fresh size-adaptive review candidates: '+out
