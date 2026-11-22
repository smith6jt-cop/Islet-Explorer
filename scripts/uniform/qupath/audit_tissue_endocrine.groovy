// Export the original cell and endocrine outlines by UUID without editing the source project.
import qupath.lib.io.GsonTools
import qupath.lib.projects.ProjectIO
import java.awt.image.BufferedImage

def argv = binding.getVariable('args') as List
def project = ProjectIO.loadProject(new File(argv[0], 'project.qpproj'), BufferedImage.class)
def folder = new File(argv[1])
def gson = GsonTools.getInstance(false)
def classes = ['Islet', 'EndocrineCluster', 'EndocrineSingle'] as Set
project.getImageList().each { entry ->
    def hierarchy = entry.readHierarchy()
    def output = new File(folder, entry.getImageName() + '.jsonl')
    output.withWriter('UTF-8') { writer ->
        def emit = { obj, kind ->
            writer.writeLine(gson.toJson([kind:kind, id:obj.getID().toString(),
                name:obj.getPathClass()?.toString(), parent:obj.getParent()?.getPathClass()?.toString(),
                x:obj.getROI().getCentroidX(), y:obj.getROI().getCentroidY(),
                geometry:obj.getROI().getGeometry()]))
        }
        hierarchy.getAnnotationObjects().findAll { classes.contains(it.getPathClass()?.toString()) }.each { emit(it,'annotation') }
        hierarchy.getDetectionObjects().findAll { it.getPathClass()?.toString() == 'Tissue cell - endocrine marker only' }.each { emit(it,'cell') }
    }
    println 'Read ' + entry.getImageName()
    System.gc()
}
