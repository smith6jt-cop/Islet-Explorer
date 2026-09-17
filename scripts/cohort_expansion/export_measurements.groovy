import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.io.GsonTools
import qupath.lib.gui.measure.ObservableMeasurementTableData
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.security.MessageDigest
import java.util.zip.GZIPOutputStream

def output = new File(System.getenv('ISLET_EXPORT_DIR'))
String caseId = System.getenv('ISLET_CASE_ID')
def cells = getCellObjects()
def annotations = getAnnotationObjects()
if (cells.isEmpty()) throw new IllegalStateException('No tissue cells to export')
def tissue = annotations.findAll { it.getPathClass()?.toString() == 'Tissue' }
if (tissue.isEmpty()) throw new IllegalStateException('No Tissue region was segmented')
def panel = ['B3TUBB', 'BActin', 'BCatenin', 'CD163', 'CD20', 'CD31', 'CD34', 'CD35', 'CD3e', 'CD4',
             'CD44', 'CD45', 'CD68', 'CD8a', 'CHGA', 'CK19', 'ColIV', 'ECAD', 'GAP43', 'GCG', 'HLADR',
             'INS', 'KRT14', 'LGALS3', 'NaKATPase', 'PDPN', 'PGP9.5', 'SMA', 'SST', 'VIM', 'panCK']
def required = panel.collect { "Cell: ${it}: Mean".toString() }
def names = cells[0].getMeasurementList().getNames()
def missing = required.findAll { !names.contains(it) }
if (!missing.isEmpty()) throw new IllegalStateException('Missing whole-cell marker measurements: ' + missing)
def calibration = getCurrentServer().getPixelCalibration()
ObjectMeasurements.addShapeMeasurements(annotations, calibration)
// The default measurement exporter rounds values according to their magnitude.
// In particular, rounded centroids can move cells across an islet boundary.
// Preserve each numeric value with Java's round-trip double representation.
def quoteField = { value ->
    String text = value == null ? '' : value.toString()
    (text.contains('\t') || text.contains('\n') || text.contains('\r') || text.contains('"')) ?
        '"' + text.replace('"', '""') + '"' : text
}
def writeMeasurements = { objects, File path ->
    def table = new ObservableMeasurementTableData()
    table.setImageData(getCurrentImageData(), objects)
    def columns = table.getAllNames().toList()
    def numeric = columns.collect { table.isNumericMeasurement(it) }
    path.withWriter('UTF-8') { writer ->
        writer.writeLine(columns.collect { quoteField(it) }.join('\t'))
        objects.each { object ->
            def values = columns.withIndex().collect { column, index ->
                numeric[index] ? Double.toString(table.getNumericValue(object, column)) : table.getStringValue(object, column)
            }
            writer.writeLine(values.collect { quoteField(it) }.join('\t'))
        }
    }
}
writeMeasurements(cells, new File(output, 'cells.tsv'))
writeMeasurements(annotations, new File(output, 'annotations.tsv'))
def geojson = [type:'FeatureCollection', features:annotations]
new File(output, 'annotations.geojson.gz').withOutputStream { stream ->
    new GZIPOutputStream(stream).withWriter('UTF-8') { writer ->
        GsonTools.getInstance().toJson(geojson, writer)
    }
}
def safeTag = getProjectEntry().getImageName().replaceAll('[^A-Za-z0-9._-]', '_')
for (kind in ['composition', 'markers', 'targets']) {
    def source = new File(buildFilePath(PROJECT_BASE_DIR, 'results', "islet_${kind}_${safeTag}.tsv"))
    if (!source.isFile()) throw new IllegalStateException('Missing regional export: ' + source)
    Files.copy(source.toPath(), new File(output, "islet_${kind}_${caseId}.tsv").toPath(), StandardCopyOption.REPLACE_EXISTING)
}
def sha256 = { File file ->
    def digest = MessageDigest.getInstance('SHA-256')
    file.withInputStream { stream ->
        byte[] buffer = new byte[1024 * 1024]
        int count
        while ((count = stream.read(buffer)) > 0) digest.update(buffer, 0, count)
    }
    digest.digest().encodeHex().toString()
}
long written = -1
new File(output, 'cells.tsv').eachLine { written++ }
if (written != cells.size()) throw new IllegalStateException("Cell export row count ${written} differs from hierarchy ${cells.size()}")
def manifest = [case_id:Integer.parseInt(caseId).toString(), source_filename:getProjectEntry().getImageName(),
                cell_count:cells.size(), annotation_count:annotations.size(), tissue_region_count:tissue.size(),
                segmentation_region:'all Tissue annotations', pixel_size_um_x:calibration.getPixelWidthMicrons(),
                pixel_size_um_y:calibration.getPixelHeightMicrons(), coordinate_units:'micrometres',
                cells_sha256:sha256(new File(output, 'cells.tsv')), annotations_sha256:sha256(new File(output, 'annotations.tsv')),
                marker_means:required, numeric_encoding:'java_double_to_string',
                exported_at:java.time.Instant.now().toString()]
new File(output, 'export_manifest.json').text = GsonTools.getInstance(true).toJson(manifest)
println "Exported ${cells.size()} tissue cells and ${annotations.size()} annotations for ${caseId}"
