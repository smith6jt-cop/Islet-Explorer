/*
 * export_tissue.groovy -- CODEX uniform pipeline (Islet-Explorer, 2026-10-15). Read-only: never saves.
 *
 * Runs against a SOURCE project (the original IsletExplorer project for the original 15 donors, the six-donor
 * review project for the additions). Exports the image's existing Tissue annotation(s) and any exclusion
 * annotations (classes starting Exclude or Ignore, or Artifact) unchanged, as a GeoJSON FeatureCollection in
 * the image's own pixel coordinates, for import_tissue.groovy (owner decision 2026-10-15: reuse the existing
 * Tissue outlines). Nothing else is exported: no cells, islets or structures.
 *
 *   args[0] = output .geojson path (written as <path>.part, re-read and checked, then renamed)
 */
import qupath.lib.io.PathIO

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: export_tissue.groovy <output.geojson>')
def outFile = new File(args[0])
// PathIO picks the reader from the extension, so the temporary file must still end in .geojson
def partFile = new File(outFile.getPath().replaceAll(/\.geojson$/, '') + '.part.geojson')
outFile.getParentFile().mkdirs()
partFile.delete()

def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def isExclusion = { String c -> c.startsWith('Exclude') || c.startsWith('Ignore') || c == 'Artifact' }
def annotations = getAnnotationObjects()
def tissue = annotations.findAll { clsOf(it) == 'Tissue' && it.getROI()?.isArea() }
def exclusions = annotations.findAll { isExclusion(clsOf(it)) && it.getROI()?.isArea() }
def cal = getCurrentServer().getPixelCalibration()
double um2 = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()

def details = [image: getProjectEntry()?.getImageName(), output: outFile.getPath(),
               tissue: tissue.size(), tissue_area_mm2: (tissue.sum { it.getROI().getArea() * um2 } ?: 0d) / 1e6,
               exclusions: exclusions.countBy { clsOf(it) },
               all_annotation_classes: annotations.countBy { clsOf(it) },
               pixel_um: [cal.getPixelWidthMicrons(), cal.getPixelHeightMicrons()]]
boolean ok = false
String message
if (tissue.isEmpty()) {
    message = 'no Tissue annotation in the source image'
} else {
    exportObjectsToGeoJson(tissue + exclusions, partFile.getPath(), 'FEATURE_COLLECTION')
    def back = PathIO.readObjects(partFile.toPath())
    int backTissue = back.count { clsOf(it) == 'Tissue' }
    double backArea = back.findAll { clsOf(it) == 'Tissue' }.sum { it.getROI().getArea() * um2 } ?: 0d
    ok = back.size() == tissue.size() + exclusions.size() && backTissue == tissue.size() &&
            Math.abs(backArea / 1e6 - details.tissue_area_mm2) < 1e-6
    message = ok ? "${tissue.size()} Tissue (${String.format('%.3f', details.tissue_area_mm2)} mm2) + ${exclusions.size()} exclusion object(s) exported" :
            "re-read mismatch: ${back.size()} objects, ${backTissue} Tissue, ${backArea / 1e6} mm2"
    if (ok && !partFile.renameTo(outFile)) { ok = false; message = "could not rename ${partFile} to ${outFile}" }
}
println "[export_tissue] ${message}"
V3_RESULT = [ok: ok, changed: false, message: message.toString(), details: details]
