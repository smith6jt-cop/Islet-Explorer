// Run an existing endocrine-only pixel classifier on raw pixels; never save source image data.
// Args: input specification JSON, classifier JSON, fresh output directory.
import qupath.lib.images.ImageData
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder
import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.images.servers.ImageChannel
import qupath.lib.classifiers.pixel.PixelClassificationImageServer
import qupath.lib.regions.RegionRequest
import qupath.lib.regions.ImagePlane
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.opencv.ml.pixel.PixelClassifierTools
import org.locationtech.jts.geom.Geometry
import java.awt.image.BufferedImage
import java.security.MessageDigest
import java.nio.file.Files
import javax.imageio.ImageIO

def argv=this.binding.hasVariable('args') ? this.binding.getVariable('args') as List : []
if(argv.size()!=3) throw new IllegalArgumentException('Provide specification, classifier and fresh output directory')
def specFile=new File(argv[0].toString()).canonicalFile
def classifierFile=new File(argv[1].toString()).canonicalFile
def out=new File(argv[2].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve an existing pilot directory')
def gson=GsonTools.getInstance(true)
def spec=gson.fromJson(specFile.text,Map.class)
def classifierDefinition=gson.fromJson(classifierFile.text,Map.class)
def classifierChannels=classifierDefinition.op.colorTransforms.collect{it.channelName.toString()}
if(classifierChannels.toSet()!=['INS','GCG','SST'].toSet() || classifierChannels.size()!=3)
    throw new IllegalStateException('This pilot requires the endocrine-only INS/GCG/SST classifier')
def sha={File file ->
    def digest=MessageDigest.getInstance('SHA-256')
    file.withInputStream{stream ->
        byte[] bytes=new byte[1048576];int n
        while((n=stream.read(bytes))>0) digest.update(bytes,0,n)
    }
    digest.digest().encodeHex().toString()
}
def source=new File(spec.source_image.toString()).canonicalFile
long sourceBytes=source.length(),sourceModified=source.lastModified()
def server=new BioFormatsServerBuilder().buildServer(source.toURI())
if(server.getWidth()!=spec.image_shape_cyx[2].intValue() || server.getHeight()!=spec.image_shape_cyx[1].intValue() ||
    server.nChannels()!=spec.image_shape_cyx[0].intValue()) throw new IllegalStateException('Source dimensions differ')
def originalChannels=server.getMetadata().getChannels()
def channels=spec.source_channel_names.withIndex().collect{name,i ->ImageChannel.getInstance(name.toString(),originalChannels[i].getColor())}
server.setMetadata(new ImageServerMetadata.Builder(server.getMetadata()).channels(channels)
    .pixelSizeMicrons(spec.pixel_size_um.doubleValue(),spec.pixel_size_um.doubleValue()).build())
def imageData=new ImageData<BufferedImage>(server,ImageData.ImageType.FLUORESCENCE)
def classifier=PixelClassifiers.readClassifier(classifierFile.toPath())
def prediction=new PixelClassificationImageServer(imageData,classifier)
double downsample=prediction.getDownsampleForResolution(0)
out.mkdirs()
Files.copy(classifierFile.toPath(),new File(out,'endocrine_classifier.json').toPath())
long started=System.nanoTime()
println "Classifying the full ${spec.case_id} image with ${classifierChannels}; downsample ${downsample}"
def request=RegionRequest.createInstance(prediction,downsample)
def mask=prediction.readRegion(request)
if(!ImageIO.write(mask,'PNG',new File(out,'pixel_classes.png'))) throw new IllegalStateException('No PNG writer')
println "Classification raster saved: ${mask.getWidth()} x ${mask.getHeight()}"
def tissue=gson.fromJson(new File(specFile.parentFile,'valid_tissue.geojson').text,Geometry.class)
def roi=GeometryTools.geometryToROI(tissue,ImagePlane.getDefaultPlane())
def objects=PixelClassifierTools.createObjectsFromPixelClassifier(prediction,
    classifier.getMetadata().getClassificationLabels(),roi,
    {r ->PathObjects.createAnnotationObject(r)} as java.util.function.Function,
    0.0d,0.0d,true)
int count=0
objects.each{object ->
    object.setName(String.format('Pixel_%s_%05d',spec.case_id,++count))
    object.getMetadata().put('Status','Unreviewed endocrine pixel-classifier output, not accepted islet geometry')
}
PathIO.exportObjectsAsGeoJSON(new File(out,'raw_pixel_regions.geojson'),objects,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
if(source.length()!=sourceBytes || source.lastModified()!=sourceModified)
    throw new IllegalStateException('Source fingerprint changed during prediction')
def metadata=[status:'pixel_classifier_pilot_unreviewed',case_id:spec.case_id,
    source_image:source.toString(),source_bytes:sourceBytes,source_mtime_ms:sourceModified,
    source_specification:specFile.toString(),source_specification_sha256:sha(specFile),
    classifier_source:classifierFile.toString(),classifier_sha256:sha(classifierFile),
    classifier_channels:classifierChannels,classifier_training_provenance:'Unknown; existing project classifier reused unchanged',
    qupath_version:qupath.lib.common.GeneralTools.getVersion(),source_pixel_size_um:spec.pixel_size_um,
    image_reader:server.getServerType(),
    prediction_downsample:downsample,prediction_width:mask.getWidth(),prediction_height:mask.getHeight(),
    source_width:server.getWidth(),source_height:server.getHeight(),
    classification_labels:classifier.getMetadata().getClassificationLabels(),
    raw_pixel_regions:objects.size(),raw_region_sha256:sha(new File(out,'raw_pixel_regions.geojson')),
    pixel_raster_sha256:sha(new File(out,'pixel_classes.png')),
    runtime_seconds:(System.nanoTime()-started)/1e9,
    original_project_modified:false,native_objects_modified:false,accepted_islets:0,
    morphology_applied:false,SAM_used:false]
new File(out,'prediction_manifest.json').text=gson.toJson(metadata)
prediction.close();imageData.close()
println gson.toJson(metadata)
