// Package the locally fitted OpenCV forest for QuPath and export a feature parity patch.
// Args: model directory, target prediction directory. No project writes.
import qupath.lib.io.GsonTools
import qupath.lib.images.ImageData
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder
import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.images.servers.ImageChannel
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.classifiers.pixel.PixelClassifierMetadata
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.RegionRequest
import qupath.opencv.ops.ImageOps
import qupath.opencv.ml.OpenCVClassifiers
import qupath.opencv.ml.pixel.PixelClassifiers
import org.bytedeco.opencv.opencv_ml.RTrees
import java.awt.image.BufferedImage

def argv=this.binding.getVariable('args') as List
if(argv.size()!=2) throw new IllegalArgumentException('Provide model and prediction directories')
def folder=new File(argv[0].toString()).canonicalFile
def gson=GsonTools.getInstance(true)
def predictionFolder=new File(argv[1].toString()).canonicalFile
def prediction=gson.fromJson(new File(predictionFolder,'prediction_manifest.json').text,Map.class)
def spec=gson.fromJson(new File(prediction.source_specification.toString()).text,Map.class)
def target=new File(predictionFolder,"Endocrine_original15_RF_${prediction.case_id}.json")
if(target.exists()) throw new IllegalStateException('Preserve packaged classifier')
def info=gson.fromJson(new File(folder,'training_manifest.json').text,Map.class)
def config=info.config
def transforms=config.markers.collect{ColorTransforms.createChannelExtractor(it.toString())}
def op=ImageOps.buildImageDataOp(transforms).appendOps(ImageOps.Core.splitMerge(
    config.gaussian_sigmas_pixels.collect{ImageOps.Filters.gaussianBlur(it.doubleValue())}),
    ImageOps.Core.subtract((prediction.normalization.offset*2).collect{it.doubleValue()} as double[]),
    ImageOps.Core.divide((prediction.normalization.scale*2).collect{it.doubleValue()} as double[]))
def forest=RTrees.load(new File(folder,'endocrine_forest.json').toString())
if(!forest.isTrained() || forest.getVarCount()!=6) throw new IllegalStateException('Invalid forest')
def model=OpenCVClassifiers.wrapStatModel(forest)
def calibration=new PixelCalibration.Builder().pixelSizeMicrons(.7548d,.7548d).build()
def metadata=new PixelClassifierMetadata.Builder().inputResolution(calibration).inputShape(512,512)
    .classificationLabels([(0):PathClass.fromString('Ignore*'),(1):PathClass.fromString('Islet')]).build()
def classifier=PixelClassifiers.createClassifier(model,op,metadata,false)
PixelClassifiers.writeClassifier(classifier,target.toPath())
// Force deserialization before declaring this a native QuPath classifier.
def reloaded=PixelClassifiers.readClassifier(target.toPath())
def server=new BioFormatsServerBuilder().buildServer(new File(spec.source_image.toString()).toURI())
def channels=spec.source_channel_names.withIndex().collect{name,i ->ImageChannel.getInstance(name.toString(),server.getMetadata().getChannels()[i].getColor())}
server.setMetadata(new ImageServerMetadata.Builder(server.getMetadata()).channels(channels).pixelSizeMicrons(.3774d,.3774d).build())
def imageData=new ImageData<BufferedImage>(server,ImageData.ImageType.FLUORESCENCE)
def bounds=[8960,5120,512,512]
def request=RegionRequest.createInstance(server.getPath(),2d,bounds[0],bounds[1],bounds[2],bounds[3])
def mat=op.apply(imageData,request)
if(!mat.isContinuous()) throw new IllegalStateException('Non-contiguous parity patch')
def buffer=mat.data().capacity(mat.total()*mat.elemSize()).asByteBuffer()
new FileOutputStream(new File(predictionFolder,'qupath_parity_features.f32')).withCloseable { stream ->
    while(buffer.hasRemaining()) stream.channel.write(buffer)
}
new File(predictionFolder,'qupath_parity_patch.json').text=gson.toJson([
    source_specification:prediction.source_specification,
    bounds_source_px:bounds,downsample:2,shape_yxc:[mat.rows(),mat.cols(),mat.channels()],
    byte_order:java.nio.ByteOrder.nativeOrder().toString(),qupath_version:qupath.lib.common.GeneralTools.getVersion(),
    feature_names:info.feature_names,classifier_roundtrip_loaded:true])
mat.close();imageData.close()
println 'Packaged native QuPath forest; feature parity must be checked before portable classifier delivery.'
