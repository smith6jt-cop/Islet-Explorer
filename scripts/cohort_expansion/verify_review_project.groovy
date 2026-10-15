// Read native review hierarchies without opening images or saving project data.
// QuPath script --args /absolute/review/root verify_review_project.groovy
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools

def argv = this.binding.hasVariable('args') ? this.binding.getVariable('args') as List : []
if (argv.size() != 1) throw new IllegalArgumentException('Provide the review project directory')
def root = new File(argv[0].toString())
def manifest = GsonTools.getInstance().fromJson(new File(root, 'review_manifest.json').text, Map.class)
def evidence = new File(root, 'verification')
evidence.mkdirs()
def required = ['B3TUBB', 'BActin', 'BCatenin', 'CD163', 'CD20', 'CD31', 'CD34', 'CD35',
                'CD3e', 'CD4', 'CD44', 'CD45', 'CD68', 'CD8a', 'CHGA', 'CK19', 'ColIV',
                'ECAD', 'GAP43', 'GCG', 'HLADR', 'INS', 'KRT14', 'LGALS3', 'NaKATPase',
                'PDPN', 'PGP9.5', 'SMA', 'SST', 'VIM', 'panCK'].collect { "Cell: ${it}: Mean".toString() }
def results = []
manifest.donors.each { donor ->
    println "Reading native object hierarchy: ${donor.case_id}"
    def hierarchy = PathIO.readHierarchy(new File(donor.review_data.toString()))
    def cells = hierarchy.getCellObjects()
    if (cells.size() != donor.cell_count.intValue())
        throw new IllegalStateException("Native cell count mismatch: ${donor.case_id}")
    def classes = hierarchy.getAnnotationObjects().countBy { it.getPathClass()?.toString() ?: 'Unclassified' }
    if ((classes['Islet'] ?: 0) != donor.islet_count.intValue())
        throw new IllegalStateException("Native islet count mismatch: ${donor.case_id}")
    def ids = new HashSet<String>()
    int nuclei = 0, fullyMeasured = 0, unmeasured = 0, partial = 0
    new File(evidence, "${donor.case_id}-cell-ids.txt").withWriter('UTF-8') { writer ->
        cells.each { cell ->
            String id = cell.getID().toString()
            if (!ids.add(id)) throw new IllegalStateException("Duplicate cell UUID: ${donor.case_id}")
            if (cell.getROI() == null) throw new IllegalStateException('Missing cell boundary')
            if (cell.getNucleusROI() != null) nuclei++
            int finite = required.count { Double.isFinite(cell.getMeasurementList().getOrDefault(it, Double.NaN)) }
            if (finite == required.size()) fullyMeasured++
            else if (finite == 0) unmeasured++
            else partial++
            writer.writeLine(id)
        }
    }
    if (partial != 0) throw new IllegalStateException("Partially missing markers: ${donor.case_id}")
    def result = [case_id:donor.case_id, cells:cells.size(), cells_with_nucleus_roi:nuclei,
                  fully_measured_cells:fullyMeasured, wholly_unmeasured_cells:unmeasured,
                  unique_cell_ids:ids.size(), annotation_classes:classes,
                  required_marker_means:required.size(), read_native_qpdata:true]
    new File(evidence, "${donor.case_id}-native.json").text = GsonTools.getInstance(true).toJson(result)
    results.add(result)
    println "Verified ${donor.case_id}: ${cells.size()} cells, ${nuclei} nuclei, ${classes['Islet']} islets"
    // Release each hierarchy before reading the next large donor file.
    cells = null; hierarchy = null; ids = null
    System.gc()
}
new File(evidence, 'native-summary.json').text = GsonTools.getInstance(true).toJson(
    [status:'passed', checked_at_utc:java.time.Instant.now().toString(), donors:results,
     project_data_modified:false, image_pixels_loaded:false])
