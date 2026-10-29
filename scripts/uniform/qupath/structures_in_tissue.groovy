/* Redraw established structure classifiers inside current Tissue. Never saves image data.
 * args: output directory, comma-separated classes (Nerve,Capillary,Lymphatic), expected Tissue GeoJSON.
 * The owner retains edge-touching structures (2026-10-29).
 */
import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.io.PathIO

if (args == null || args.length < 3) throw new IllegalArgumentException('output directory, classes and expected Tissue required')
def output = new File(args[0]); output.mkdirs()
def requested = args[1].split(',').toList()
def names = ['Nerve':'Nerve', 'Capillary':'Capillary', 'Lymphatic':'Lymphatics']
def minimum = ['Nerve':2.0d, 'Capillary':4.0d, 'Lymphatic':4.0d]
def clsOf = { o -> o.getPathClass()?.toString() ?: '' }
def tissues = getAnnotationObjects().findAll { clsOf(it) == 'Tissue' }
if (tissues.size() != 1) throw new IllegalStateException('Expected one Tissue annotation')
def tissue = tissues[0]
def region = tissue.getROI().getGeometry()
def expectedTissues = PathIO.readObjects(new File(args[2]).toPath()).findAll { clsOf(it) == 'Tissue' }
if (expectedTissues.size() != 1)
    throw new IllegalStateException('Expected Tissue GeoJSON must contain exactly one Tissue annotation')
def expectedRegion = expectedTissues[0].getROI().getGeometry()
// A matching provenance filename is not enough: inference uses the Tissue embedded in qpdata.
// Symmetric difference allows representation/order changes while rejecting any positive-area mismatch.
if (!region.isValid() || !expectedRegion.isValid() || region.symDifference(expectedRegion).getArea() > 1e-6)
    throw new IllegalStateException('Embedded QuPath Tissue differs from the supplied reviewed outline; install it first')
for (cls in requested) {
    if (!names.containsKey(cls)) throw new IllegalArgumentException('Unknown class ' + cls)
    def previous = getAnnotationObjects().findAll { clsOf(it) == cls }
    if (!previous.isEmpty()) removeObjects(previous, true)
    selectObjects(tissues)
    createAnnotationsFromPixelClassifier(names[cls], minimum[cls], -1.0, 'SPLIT')
    def drawn = getAnnotationObjects().findAll { clsOf(it) == cls }
    // Exact intersection prevents raster conversion from leaking outside the final Tissue ROI.
    def clipped = drawn.collect { o ->
        def g = o.getROI().getGeometry().intersection(region)
        if (g.isEmpty()) return null
        return PathObjects.createAnnotationObject(GeometryTools.geometryToROI(g, o.getROI().getImagePlane()), o.getPathClass())
    }.findAll { it != null }
    exportObjectsToGeoJson(clipped, new File(output, cls + '.geojson').getPath(), 'FEATURE_COLLECTION')
    println "[structures] ${cls}: ${clipped.size()} objects; edge-touching objects retained"
}
