/*
 * lymphatics_in_tissue.groovy -- CODEX uniform pipeline (Islet-Explorer, 2026-10-27). Read-only: never saves.
 *
 * The additions' Lymphatic objects are the original project's `Lymphatics` pixel classifier's own objects
 * (typing/structure_thresholds.py KEEP_CLASSIFIER). They were drawn inside the review project's Tissue, which the
 * owner found defective; this draws them again inside the image's current (redrawn) Tissue, as the original project's
 * 4_Islet_Cells.groovy drew them: createAnnotationsFromPixelClassifier("Lymphatics", 4.0, -1.0, "SPLIT") with Tissue
 * selected, then every Lymphatic object touching the 10 um ring at the image edge removed.
 *
 *   args[0] = output .geojson (FeatureCollection of the Lymphatic objects, full-resolution pixel coordinates)
 */
import qupath.lib.roi.ROIs

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: lymphatics_in_tissue.groovy <out.geojson>')
def outFile = new File(args[0])
outFile.getParentFile().mkdirs()

def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def tissue = getAnnotationObjects().findAll { clsOf(it) == 'Tissue' }
if (tissue.size() != 1) throw new IllegalStateException("[lymphatics] expected one Tissue annotation, found ${tissue.size()}")
def before = getAnnotationObjects().findAll { clsOf(it) == 'Lymphatic' }
if (!before.isEmpty()) removeObjects(before, true)

selectObjectsByClassification('Tissue')
createAnnotationsFromPixelClassifier('Lymphatics', 4.0, -1.0, 'SPLIT')

def server = getCurrentServer()
def cal = server.getPixelCalibration()
double mx = 10.0 / cal.getPixelWidthMicrons(), my = 10.0 / cal.getPixelHeightMicrons()
def drawn = getAnnotationObjects().findAll { clsOf(it) == 'Lymphatic' }
def inner = ROIs.createRectangleROI(mx, my, server.getWidth() - 2 * mx, server.getHeight() - 2 * my, null).getGeometry()
def kept = drawn.findAll { inner.contains(it.getROI().getGeometry()) }
exportObjectsToGeoJson(kept, outFile.getPath(), 'FEATURE_COLLECTION')
double um2 = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()
double area = (kept.sum { it.getROI().getArea() } ?: 0d) * um2
double tissueArea = tissue[0].getROI().getArea() * um2
println "[lymphatics] drawn ${drawn.size()}, kept ${kept.size()} (edge ring removed ${drawn.size() - kept.size()}); " +
        "area ${String.format('%.0f', area)} um2 = ${String.format('%.3f', 100 * area / tissueArea)} % of Tissue -> ${outFile}"
