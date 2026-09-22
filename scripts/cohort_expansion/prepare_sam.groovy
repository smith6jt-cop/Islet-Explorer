import qupath.lib.display.ImageDisplay
import qupath.lib.gui.images.servers.RenderedImageServer
import qupath.lib.io.GsonTools
import qupath.lib.regions.RegionRequest
import javax.imageio.ImageIO

def server = getCurrentServer()
def output = new File(System.getenv('ISLET_EXPORT_DIR'))
output.mkdirs()
def islets = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
if (islets.isEmpty()) throw new IllegalStateException('No candidate islets from pixel classifier')
def display = ImageDisplay.create(getCurrentImageData())
def wanted = ['SST', 'GCG', 'INS'] as Set
display.availableChannels().each { channel ->
    String name = channel.getName().replaceAll('\\s*\\(C\\d+\\)\\s*$', '')
    display.setChannelSelected(channel, wanted.contains(name))
}
def rendered = new RenderedImageServer.Builder(getCurrentImageData()).renderer(display).build()
double downsample = 6.0
def request = RegionRequest.createInstance(rendered.getPath(), downsample, 0, 0, server.getWidth(), server.getHeight(), 0, 0)
def image = rendered.readRegion(request)
if (!ImageIO.write(image, 'PNG', new File(output, 'sam_render.png'))) throw new IOException('PNG writer unavailable')
def boxes = []
islets.eachWithIndex { object, index ->
    def roi = object.getROI()
    double x = Math.max(0d, roi.getBoundsX() - 16d)
    double y = Math.max(0d, roi.getBoundsY() - 16d)
    double width = Math.min(server.getWidth() - x, roi.getBoundsWidth() + 32d)
    double height = Math.min(server.getHeight() - y, roi.getBoundsHeight() + 32d)
    if (width >= 4 && height >= 4) {
        def boxROI = qupath.lib.roi.ROIs.createRectangleROI(x, y, width, height, roi.getImagePlane())
        def region = RegionRequest.createInstance(rendered.getPath(), downsample, boxROI)
        boxes << [index:index, source_object_id:object.getID().toString(),
                  bbox:[(int)(region.getMinX()/downsample), (int)(region.getMinY()/downsample),
                        (int)Math.round(region.getMaxX()/downsample), (int)Math.round(region.getMaxY()/downsample)]]
    }
}
if (boxes.isEmpty()) throw new IllegalStateException('No valid SAM prompt boxes')
def metadata = [source_filename:getProjectEntry().getImageName(), downsample:downsample,
                width:image.getWidth(), height:image.getHeight(), model:'vit_b', multimask_output:false,
                display:GsonTools.getInstance().fromJson(display.toJSON(false), Object.class), boxes:boxes]
new File(output, 'sam_prompts.json').text = GsonTools.getInstance(true).toJson(metadata)
println "Prepared ${boxes.size()} SAM prompts on ${image.getWidth()}x${image.getHeight()} rendered pixels"
rendered.close()
