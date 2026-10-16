/*
 * import_tissue.groovy -- CODEX uniform pipeline (Islet-Explorer, 2026-10-15). First step on a fresh entry of the
 * uniform project; the runner's footer saves the image only if V3_RESULT.ok.
 *
 * Imports the Tissue (and exclusion) annotations that export_tissue.groovy wrote from the source project, with
 * their outlines unchanged (owner decision 2026-10-15: reuse the existing Tissue outlines). Then:
 *   - exclusion classes (Exclude*, Ignore*) become 'Artifact', the class Senior's scripts 3, 4 and 8 handle;
 *   - several Tissue pieces are merged into ONE Tissue annotation with the same outline, because Senior's cell
 *     detection (0_cell_detection.groovy) and preflight require exactly one;
 *   - image type FLUORESCENCE and pixel size 0.3774 um, as every later step sets.
 * Refuses an entry that already holds detections: re-importing would orphan them.
 *
 *   args[0] = the exported .geojson
 */
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.roi.RoiTools

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: import_tissue.groovy <tissue.geojson>')
def inFile = new File(args[0])
if (!inFile.isFile()) throw new IllegalArgumentException("[import_tissue] not found: ${inFile}")

setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.3774, 0.3774)
def hierarchy = getCurrentHierarchy()
if (!hierarchy.getDetectionObjects().isEmpty())
    throw new IllegalStateException("[import_tissue] ${hierarchy.getDetectionObjects().size()} detections exist: import runs on a fresh entry only")
def previous = new ArrayList(hierarchy.getAnnotationObjects())
if (!previous.isEmpty()) {
    println "[import_tissue] replacing ${previous.size()} annotation(s) from a previous import"
    removeObjects(previous, false)
}

def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def objects = PathIO.readObjects(inFile.toPath())
def tissuePieces = objects.findAll { clsOf(it) == 'Tissue' }
def exclusions = objects.findAll { clsOf(it) != 'Tissue' }
def renamed = exclusions.countBy { clsOf(it) }
exclusions.each { it.setPathClass(getPathClass('Artifact')) }

def cal = getCurrentServer().getPixelCalibration()
double um2 = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()
double areaBefore = (tissuePieces.sum { it.getROI().getArea() } ?: 0d) as double
def tissue = tissuePieces
boolean merged = false
if (tissuePieces.size() > 1) {
    def roi = RoiTools.union(tissuePieces.collect { it.getROI() })
    tissue = [PathObjects.createAnnotationObject(roi, getPathClass('Tissue'))]
    merged = true
}
addObjects(tissue + exclusions)
resolveHierarchy()

def tissueNow = getAnnotationObjects().findAll { clsOf(it) == 'Tissue' }
double areaNow = (tissueNow.sum { it.getROI().getArea() } ?: 0d) as double
// a merge of overlapping pieces may shrink the summed area; it can never grow it
boolean ok = tissueNow.size() == 1 && areaNow > 0 && areaNow <= areaBefore * (1 + 1e-9)
def details = [source: inFile.getPath(), tissue_pieces_in: tissuePieces.size(), merged: merged,
               tissue_area_mm2: areaNow * um2 / 1e6, pieces_area_mm2: areaBefore * um2 / 1e6,
               exclusions_as_artifact: renamed, pixel_um: 0.3774]
String message = ok ? "1 Tissue (${String.format('%.3f', areaNow * um2 / 1e6)} mm2${merged ? ', merged from ' + tissuePieces.size() + ' pieces' : ''}) + ${exclusions.size()} Artifact" :
        "expected exactly one Tissue after import, found ${tissueNow.size()}"
println "[import_tissue] ${message}"
V3_RESULT = [ok: ok, changed: true, message: message.toString(), details: details]
