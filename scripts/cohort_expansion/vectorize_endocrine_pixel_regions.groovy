// Convert externally classified integer-grid rasters using QuPath's native tracing.
// The uniform boundary raster is a historical comparator. For current delivery,
// follow this step with build_size_adaptive_islets.groovy on the raw foreground.
// Args: prediction directory. Never load or save a source project's image data.
import qupath.lib.analysis.images.ContourTracing
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.regions.RegionRequest
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import org.locationtech.jts.geom.Geometry
import java.security.MessageDigest
import javax.imageio.ImageIO

def argv=this.binding.getVariable('args') as List
if(argv.size()!=1) throw new IllegalArgumentException('Provide prediction directory')
def folder=new File(argv[0].toString()).canonicalFile
def gson=GsonTools.getInstance(true)
def manifestFile=new File(folder,'prediction_manifest.json')
def manifest=gson.fromJson(manifestFile.text,Map.class)
if(new File(folder,'vectorization_manifest.json').exists()) throw new IllegalStateException('Preserve existing vectorization')
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def tissueFile=new File(new File(manifest.source_specification.toString()).parentFile,'valid_tissue.geojson')
def tissue=gson.fromJson(tissueFile.text,Geometry.class)
def request=RegionRequest.createInstance('endocrine-pixel-refinement',manifest.prediction_downsample.doubleValue(),
    0,0,manifest.source_width.intValue(),manifest.source_height.intValue())
def outputs=[]
['pixel_classes','boundary_regions'].each{kind ->
    def rasterFile=new File(folder,kind+'.png')
    if(sha(rasterFile)!=manifest.raster_sha256[kind+'.png']) throw new IllegalStateException('Prediction raster changed')
    def image=ImageIO.read(rasterFile)
    if(image.width!=manifest.prediction_width.intValue() || image.height!=manifest.prediction_height.intValue())
        throw new IllegalStateException('Raster dimensions differ')
    def geometry=ContourTracing.createTracedGeometry(image.getRaster(),1d,1d,0,request)
    if(geometry==null) throw new IllegalStateException('No foreground geometry')
    double before=geometry.area
    def clipped=geometry.intersection(tissue)
    def objects=[]
    for(int i=0;i<clipped.numGeometries;i++) {
        def piece=clipped.getGeometryN(i)
        if(piece.isEmpty() || piece.area==0) continue
        def roi=GeometryTools.geometryToROI(piece,ImagePlane.getDefaultPlane())
        def object=PathObjects.createAnnotationObject(roi,PathClass.fromString(
            kind=='boundary_regions' ? 'Endocrine pixel candidate' : 'Endocrine pixel foreground'))
        object.setName(String.format('%s_%s_%05d',kind=='boundary_regions' ? 'PixelIslet':'Pixel',manifest.case_id,objects.size()+1))
        object.setColor(0,255,255)
        object.getMeasurementList().put('Pixel candidate: area um2',piece.area*Math.pow(manifest.source_pixel_size_um.doubleValue(),2))
        object.getMetadata().put('Review status','Unreviewed reference-trained endocrine pixel region; not accepted islet geometry')
        object.getMetadata().put('Training donors','Original 15 donors; HPAP078 and all six additions excluded')
        object.getMetadata().put('Segmentation method','INS/GCG/SST random forest, 2 um closing, enclosed holes up to 80 um2; no SAM prediction')
        object.getMetadata().put('Model SHA256',manifest.model_sha256.toString())
        objects.add(object)
    }
    def target=new File(folder,kind=='boundary_regions' ? 'review_candidates.geojson' : 'raw_pixel_regions.geojson')
    if(target.exists()) throw new IllegalStateException('Preserve existing vectors')
    PathIO.exportObjectsAsGeoJSON(target,objects,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
    outputs.add([kind:kind,objects:objects.size(),unclipped_area_source_px2:before,clipped_area_source_px2:clipped.area,
        geometry_valid:clipped.isValid(),path:target.canonicalPath,sha256:sha(target)])
    println "${kind}: ${objects.size()} native regions"
}
new File(folder,'vectorization_manifest.json').text=gson.toJson([
    status:'native_raster_vectorization_complete',qupath_version:qupath.lib.common.GeneralTools.getVersion(),
    source_prediction_manifest_sha256:sha(manifestFile),tissue_sha256:sha(tissueFile),
    prediction_downsample:manifest.prediction_downsample,outputs:outputs,
    original_project_modified:false,source_native_objects_modified:false,accepted_islets:0])
