/* Sanitize a fresh, private copy for native image review. Never run against the source project.
 * args[0]: private JSON assignments (object_id, sample_id, blind_uuid); args[1]: private report path.
 */
import com.google.gson.Gson
import qupath.lib.objects.classes.PathClass
import org.locationtech.jts.io.WKBWriter

if (args == null || args.length != 2) throw new IllegalArgumentException('assignment JSON and report required')
def rows = new Gson().fromJson(new File(args[0]).getText('UTF-8'), List)
def assignments = rows.collectEntries { [(it.object_id): it] }
if (assignments.size() != rows.size()) throw new IllegalStateException('Duplicate source IDs in blind assignments')
if (rows.collect { it.sample_id }.toSet().size() != rows.size() || rows.collect { it.blind_uuid }.toSet().size() != rows.size())
    throw new IllegalStateException('Duplicate blind IDs or names')
def hierarchy = getCurrentHierarchy()
def cells = hierarchy.getDetectionObjects().findAll { assignments.containsKey(it.getID().toString()) }
if (cells.size() != rows.size()) throw new IllegalStateException('Some sampled cells are missing from native project')
if (cells.any { !it.isCell() }) throw new IllegalStateException('Sample contains a non-cell detection')
def tissues = hierarchy.getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Tissue' }
if (tissues.size() != 1) throw new IllegalStateException('Expected one Tissue annotation')
def tissue = tissues[0]
def geometry = { obj ->
    [new WKBWriter().write(obj.getROI().getGeometry()).encodeHex().toString(),
     obj.isCell() && obj.getNucleusROI() != null ?
         new WKBWriter().write(obj.getNucleusROI().getGeometry()).encodeHex().toString() : null]
}
def before = cells.collectEntries { [(assignments[it.getID().toString()].sample_id): geometry(it)] }
def tissueGeometry = geometry(tissue)
def neutral = PathClass.fromString('Review cell', 0x808080)
cells.each { cell ->
    def row = assignments[cell.getID().toString()]
    cell.getParent()?.removeChildObject(cell)
    cell.clearChildObjects()
    cell.setID(java.util.UUID.fromString(row.blind_uuid))
    cell.setName(row.sample_id)
    cell.setPathClass(neutral, Double.NaN)
    cell.setColor(null)
    cell.getMetadata().clear()
    cell.getMeasurementList().clear()
}
tissue.clearChildObjects()
tissue.setName(null)
tissue.setDescription(null)
tissue.setPathClass(PathClass.fromString('Tissue', 0x669999), Double.NaN)
tissue.setColor(null)
tissue.getMetadata().clear()
tissue.getMeasurementList().clear()
hierarchy.clearAll()
hierarchy.getRootObject().getMetadata().clear()
hierarchy.getRootObject().getMeasurementList().clear()
hierarchy.getRootObject().setName(null)
hierarchy.getRootObject().setColor(null)
hierarchy.addObject(tissue)
tissue.addChildObjects(cells)
hierarchy.getSelectionModel().clearSelection()
hierarchy.fireHierarchyChangedEvent(this)
def imageData = getCurrentImageData()
new ArrayList(imageData.getProperties().keySet()).each { imageData.removeProperty(it) }
imageData.getHistoryWorkflow().clear()
if (geometry(tissue) != tissueGeometry || cells.any { geometry(it) != before[it.getName()] })
    throw new IllegalStateException('Cell/nucleus/Tissue geometry changed during blinding')
if (hierarchy.getDetectionObjects().size() != rows.size() || hierarchy.getAnnotationObjects().size() != 1)
    throw new IllegalStateException('Unexpected cells or annotations remain after blinding')
if (cells.any { !it.getMetadata().isEmpty() || !it.getMeasurementList().isEmpty() ||
                it.getPathClass().toString() != 'Review cell' } ||
    !imageData.getProperties().isEmpty() || !imageData.getHistoryWorkflow().getSteps().isEmpty())
    throw new IllegalStateException('Prediction metadata or workflow history remains')
def entry = getProjectEntry()
if (entry == null) throw new IllegalStateException('Blind import requires a native project entry')
entry.saveImageData(imageData)
def report = [passed: true, image: entry.getImageName(), samples: cells.size(), tissue_annotations: 1,
              geometry_preserved: true, metadata_cleared: true, measurements_cleared: true,
              workflow_cleared: true, accuracy_computed: false]
new File(args[1]).setText(new Gson().toJson(report), 'UTF-8')
println "[blind review] saved ${entry.getImageName()}: ${cells.size()} neutral samples"
