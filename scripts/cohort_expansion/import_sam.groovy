import qupath.lib.io.PathIO
import java.nio.file.Paths

def objects = PathIO.readObjects(Paths.get(System.getenv('ISLET_EXPORT_DIR'), 'sam_islets.geojson'))
def cleaned = objects.findAll { object ->
    def roi = object.getROI()
    roi != null && roi.getBoundsWidth() >= 2 && roi.getBoundsHeight() >= 2 && roi.getArea() >= 16
}
if (cleaned.isEmpty()) throw new IllegalStateException('No valid SAM masks after tiny-mask filtering')
cleaned.each {
    it.setPathClass(getPathClass('Islet')); it.setName('Islet')
    // Geometric validity and SAM mask confidence do not establish islet identity.
    // Retain candidates for measurement and require an explicit review before fitting.
    it.getMeasurementList().put('Islet: identity reviewed', 0.0)
}
def originals = getAnnotationObjects().findAll { it.getPathClass()?.toString() in ['Islet', 'Box'] }
removeObjects(originals, true)
addObjects(cleaned)
fireHierarchyUpdate()
println "Imported ${cleaned.size()} SAM islets; removed ${objects.size()-cleaned.size()} tiny masks"
println 'Imported regions are unreviewed islet candidates; review endocrine channels and tissue context before analytical acceptance.'
