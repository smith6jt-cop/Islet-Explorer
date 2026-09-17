import qupath.lib.display.ImageDisplay
import qupath.lib.gui.images.servers.RenderedImageServer
import qupath.lib.io.GsonTools
import qupath.lib.regions.RegionRequest
import javax.imageio.ImageIO

def output = new File(System.getenv('ISLET_EXPORT_DIR'), 'qc')
output.mkdirs()
def server = getCurrentServer()
def cells = getCellObjects()
def candidates = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() == 'Islet' && it.getMeasurementList().getOrDefault('Islet: cells total', 0) >= 20 &&
    it.getROI().getBoundsWidth() <= 600 && it.getROI().getBoundsHeight() <= 600
}.sort { it.getMeasurementList().getOrDefault('Islet: cells total', 0) }
if (candidates.isEmpty()) candidates = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }.sort { it.getROI().getArea() }
if (candidates.isEmpty() || cells.isEmpty()) throw new IllegalStateException('No segmented islet/cells for crop review')
def islet = candidates[(int)(candidates.size()/2)]
int width = Math.min(1024, server.getWidth())
int height = Math.min(1024, server.getHeight())
int x = Math.max(0, Math.min(server.getWidth()-width, (int)Math.round(islet.getROI().getCentroidX()-width/2)))
int y = Math.max(0, Math.min(server.getHeight()-height, (int)Math.round(islet.getROI().getCentroidY()-height/2)))
def crop = qupath.lib.roi.ROIs.createRectangleROI(x, y, width, height, islet.getROI().getImagePlane())
def display = ImageDisplay.create(getCurrentImageData())
def wanted = ['DAPI', 'INS', 'GCG'] as Set
display.availableChannels().each { channel ->
    display.setChannelSelected(channel, wanted.contains(channel.getName().replaceAll('\\s*\\(C\\d+\\)\\s*$', '')))
}
def rendered = new RenderedImageServer.Builder(getCurrentImageData()).renderer(display).build()
def request = RegionRequest.createInstance(rendered.getPath(), 1.0, x, y, width, height, 0, 0)
def image = rendered.readRegion(request)
ImageIO.write(image, 'PNG', new File(output, 'DAPI_INS_GCG.png'))
def cropGeometry = crop.getGeometry()
def cropCells = cells.findAll { it.getROI().getGeometry().intersects(cropGeometry) }
def cropAnnotations = getAnnotationObjects().findAll {
    it.getPathClass()?.toString() in ['Islet', 'IsletExpanded'] && it.getROI().getGeometry().intersects(cropGeometry)
}
new File(output, 'objects.geojson').text = GsonTools.getInstance().toJson([type:'FeatureCollection', features:cropCells+cropAnnotations])
def metadata = [case_id:System.getenv('ISLET_CASE_ID'), source_filename:getProjectEntry().getImageName(),
                x_px:x, y_px:y, width_px:width, height_px:height, downsample:1.0,
                pixel_size_um:server.getPixelCalibration().getAveragedPixelSizeMicrons(),
                representative_islet:islet.getName(), selection:'Median cell count among islets with at least20cells and bounds<=600pixels; medianarea fallback',
                cell_count:cropCells.size(), annotation_count:cropAnnotations.size(),
                display:GsonTools.getInstance().fromJson(display.toJSON(false), Object.class)]
new File(output, 'crop.json').text = GsonTools.getInstance(true).toJson(metadata)
println "QC crop exported around ${islet.getName()}: ${cropCells.size()} cells"
rendered.close()
