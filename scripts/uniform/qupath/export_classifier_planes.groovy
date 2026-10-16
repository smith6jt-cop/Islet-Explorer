/*
 * export_classifier_planes.groovy -- CODEX uniform pipeline (Islet-Explorer, 2026-10-15). Read-only: never saves.
 *
 * Writes the islet classifier's input channels at the classifier's resolution, as QuPath itself reads them, for
 * islet_classifier_fit.py --planes. Senior's fit read pyramid level 1 of an Akoya .qptiff directly; these CODEX
 * images come in several formats, so QuPath does the reading here: the fit then sees exactly the pixels the
 * classifier will see (downsample 2 = 0.7548 um, the classifier's resolution).
 *
 *   args[0] = output .tif (float32 ImageJ stack, one slice per channel); a .json sidecar is written beside it
 *   args[1] = downsample (default 2.0)
 */
import qupath.lib.regions.RegionRequest
import ij.ImagePlus
import ij.ImageStack
import ij.io.FileSaver
import ij.process.FloatProcessor
import qupath.lib.io.GsonTools

List<String> CHANNELS = ['CHGA', 'INS', 'GCG', 'SST']   // islet_classifier_fit.CHANNELS
int TILE_OUT = 2048                                       // tile edge in output pixels

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: export_classifier_planes.groovy <out.tif> [downsample]')
def outFile = new File(args[0])
double ds = args.length > 1 ? Double.parseDouble(args[1]) : 2.0
def partFile = new File(outFile.getPath() + '.part')
def sidecar = new File(outFile.getPath().replaceAll(/\.tiff?$/, '') + '.json')
outFile.getParentFile().mkdirs()
partFile.delete()

def server = getCurrentServer()
def names = server.getMetadata().getChannels().collect { it.getName() }
def index = CHANNELS.collect { names.indexOf(it) }
def missing = [CHANNELS, index].transpose().findAll { it[1] < 0 }.collect { it[0] }
if (!missing.isEmpty()) throw new IllegalStateException("[planes] channel(s) not in this image: ${missing}")
double pixelUm = server.getPixelCalibration().getAveragedPixelSizeMicrons()
int W = (int) Math.round(server.getWidth() / ds)
int H = (int) Math.round(server.getHeight() / ds)
float[][] planes = new float[CHANNELS.size()][W * H]
long t0 = System.currentTimeMillis()

for (int oy = 0; oy < H; oy += TILE_OUT) {
    for (int ox = 0; ox < W; ox += TILE_OUT) {
        int x = (int) Math.round(ox * ds), y = (int) Math.round(oy * ds)
        int w = (int) Math.min(Math.round(TILE_OUT * ds), server.getWidth() - x)
        int h = (int) Math.min(Math.round(TILE_OUT * ds), server.getHeight() - y)
        if (w <= 0 || h <= 0) continue
        def img = server.readRegion(RegionRequest.createInstance(server.getPath(), ds, x, y, w, h))
        def raster = img.getRaster()
        int tw = Math.min(img.getWidth(), W - ox), th = Math.min(img.getHeight(), H - oy)
        float[] buf = new float[tw * th]
        CHANNELS.eachWithIndex { c, ci ->
            raster.getSamples(0, 0, tw, th, index[ci], buf)
            for (int r = 0; r < th; r++) System.arraycopy(buf, r * tw, planes[ci], (oy + r) * W + ox, tw)
        }
    }
}

def stack = new ImageStack(W, H)
CHANNELS.eachWithIndex { c, ci -> stack.addSlice(c, new FloatProcessor(W, H, planes[ci])) }
def imp = new ImagePlus(getProjectEntry().getImageName() + ' classifier planes', stack)
boolean written = new FileSaver(imp).saveAsTiffStack(partFile.getPath())
boolean ok = written && partFile.length() > 0 && partFile.renameTo(outFile)
def meta = [channels: CHANNELS, downsample: ds, pixel_um: ds * pixelUm, image_pixel_um: pixelUm, width: W, height: H,
            image: getProjectEntry().getImageName(), server_path: server.getPath(),
            seconds: (System.currentTimeMillis() - t0) / 1000.0]
if (ok) sidecar.text = GsonTools.getInstance(true).toJson(meta)
String message = ok ? "${CHANNELS.size()} planes ${W} x ${H} at ${ds * pixelUm} um" : "write failed (${partFile})"
println "[planes] ${message}"
V3_RESULT = [ok: ok, changed: false, message: message.toString(), details: meta + [output: outFile.getPath(), sidecar: sidecar.getPath()]]
