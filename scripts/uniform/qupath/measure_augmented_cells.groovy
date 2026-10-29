/* Measure only explicit IDs written by augment_tissue.groovy; never select every detection. */
import com.google.gson.Gson
import qupath.lib.analysis.features.ObjectMeasurements
import java.util.concurrent.ForkJoinPool

def cfg = new Gson().fromJson(new File(args[0]).getText('UTF-8'), Map)
def invariant = evaluate(new File(cfg.invariants_script))
def saved = invariant.readJSON(new File(cfg.inventory))
def ids = invariant.readJSON(new File(cfg.new_ids)) as Set
def hierarchy = getCurrentHierarchy()
invariant.verify(hierarchy, saved)
if (!Collections.disjoint(ids, saved.cells.keySet())) throw new IllegalStateException('New IDs overlap original IDs')
def selected = hierarchy.getDetectionObjects().findAll { ids.contains(it.getID().toString()) }
if (selected.size() != ids.size()) throw new IllegalStateException('Some new IDs are missing from the hierarchy')
def imageData = getCurrentImageData()
def server = imageData.getServerBuilder().build()
server.setMetadata(imageData.getServer().getMetadata())
def failures = new java.util.concurrent.ConcurrentLinkedQueue()
def pool = new ForkJoinPool(24)
try {
    if (selected) {
        hierarchy.getSelectionModel().setSelectedObjects(selected, null)
        addShapeMeasurements('AREA', 'LENGTH', 'CIRCULARITY', 'SOLIDITY', 'MAX_DIAMETER', 'MIN_DIAMETER', 'NUCLEUS_CELL_RATIO')
        int tw = server.getMetadata().getPreferredTileWidth(), th = server.getMetadata().getPreferredTileHeight()
        long nx = (long) Math.ceil(server.getWidth() / (double) tw)
        selected = new ArrayList(selected)
        selected.sort { a, b ->
            def ra = a.getROI(), rb = b.getROI()
            long ta = ((long) (ra.getCentroidY() / th)) * nx + (long) (ra.getCentroidX() / tw)
            long tb = ((long) (rb.getCentroidY() / th)) * nx + (long) (rb.getCentroidX() / tw)
            ta <=> tb
        }
        pool.submit {
            selected.parallelStream().forEach { detection ->
                try {
                    ObjectMeasurements.addIntensityMeasurements(server, detection, 1.0,
                        [ObjectMeasurements.Measurements.MEAN, ObjectMeasurements.Measurements.MIN,
                         ObjectMeasurements.Measurements.MAX, ObjectMeasurements.Measurements.STD_DEV],
                        [ObjectMeasurements.Compartments.NUCLEUS, ObjectMeasurements.Compartments.CYTOPLASM,
                         ObjectMeasurements.Compartments.MEMBRANE, ObjectMeasurements.Compartments.CELL])
                } catch (Exception e) { failures.add("${detection.getID()}: ${e}") }
            }
        }.get()
    }
} finally {
    pool.shutdown(); server.close(); hierarchy.getSelectionModel().clearSelection()
}
if (failures) throw new IllegalStateException("Measurement failed for ${failures.size()} cells: ${failures.peek()}")
def missing = selected.findAll { it.getROI().getArea() * 0.3774 * 0.3774 >= 5.0 &&
    Double.isNaN(it.getMeasurementList().get('Cell: DAPI: Mean')) }
if (missing) throw new IllegalStateException("${missing.size()} new cells lack DAPI measurements")
def checks = invariant.verify(hierarchy, saved)
invariant.requireSaveWindow(cfg)
hierarchy.fireHierarchyChangedEvent(this)
V3_RESULT = [ok: true, changed: !selected.isEmpty(), message: "Measured ${selected.size()} new cells only",
             details: [measured: selected.size(), preservation: checks]]
