// Apply UUID deltas only to a separate candidate project and read back every detection.
import qupath.lib.projects.ProjectIO
import qupath.lib.objects.classes.PathClass
import qupath.lib.io.GsonTools
import java.awt.image.BufferedImage

def argv = this.binding.getVariable('args') as List
def root = new File(argv[0])
def evidence = new File(argv[1])
if (!new File(root, 'V24_REVIEW_COPY').isFile()) throw new IllegalArgumentException('Expected candidate copy')
def project = ProjectIO.loadProject(new File(root,'project.qpproj'), BufferedImage.class)
def gson = GsonTools.getInstance(false)
def catalog = project.getPathClasses().collectEntries { [(it.toString()):it] }
def usedClasses = [] as Set
def terminalReview = new File(evidence, 'terminal_classes.json')
def policy = terminalReview.isFile() ? gson.fromJson(terminalReview.getText('UTF-8'), Map) : null
def allowedDisplays = policy == null ? null : policy.detection_classes as Set
project.getImageList().each { entry ->
    def name = entry.getImageName()
    def receipt = new File(evidence,'samples/'+name+'.json')
    def lines = new File(evidence,'qupath/'+name+'.tsv').readLines('UTF-8')
    def header = lines.first().split('\t',-1)
    def changes = lines.drop(1).collectEntries { line ->
        def record = [header.toList(),line.split('\t',-1).toList()].transpose().collectEntries { [(it[0]):it[1]] }
        [(record.object_id):record]
    }
    def sampled = new File(evidence,'samples/'+name+'.txt').readLines() as Set
    def data = entry.readImageData()
    def expected = [:]
    def geometry = [:]
    int changed = 0
    data.getHierarchy().getDetectionObjects().each { cell ->
        def id = cell.getID().toString()
        def r = changes[id]
        if (r) {
            def meta = cell.getMetadata()
            if ((r.old_display && cell.getPathClass()?.toString() != r.old_display) || meta['phenotype_label'] != r.old_label) throw new IllegalStateException('Source mismatch '+id)
            def cls = catalog[r.display] ?: PathClass.getInstance(r.display)
            cell.setPathClass(cls)
            meta['phenotype_label'] = r.label
            def fields = r.containsKey('metadata_keys') ? r.metadata_keys.tokenize(';') :
                    (header.toList()-['object_id','old_label','old_display','label','display'])
            fields.each { k ->
                if (r[k]) meta[k] = r[k]
                else meta.remove(k)
            }
            changed++
        }
        if (sampled.contains(id)) geometry[id] = cell.getROI().getAllPoints().collect { [it.getX(),it.getY()] }
        expected[id] = cell.getPathClass()?.toString()
    }
    if (changed != changes.size() || geometry.size()!=sampled.size()) throw new IllegalStateException(name+' UUID mismatch')
    if (allowedDisplays != null && !allowedDisplays.containsAll(expected.values()))
        throw new IllegalStateException(name+' has unapproved detection classes '+(expected.values().toSet()-allowedDisplays))
    usedClasses.addAll(expected.values())
    usedClasses.addAll(data.getHierarchy().getAnnotationObjects().collect { it.getPathClass()?.toString() })
    entry.saveImageData(data)
    data = null
    def actual = entry.readHierarchy().getDetectionObjects()
    if (actual.size()!=expected.size() || !actual.every { obj ->
        def id=obj.getID().toString()
        expected.containsKey(id) && obj.getPathClass()?.toString()==expected[id] &&
        (!changes[id] || (obj.getMetadata()['phenotype_label']==changes[id].label &&
          (!policy || (obj.getMetadata()['excluded_reason'] ?: '') == changes[id].excluded_reason) &&
          (!changes[id].subtype_status || obj.getMetadata()['subtype_status'] == changes[id].subtype_status) &&
          (!changes[id].t_subtype_cd8_vote || obj.getMetadata()['t_subtype_cd8_vote'] == changes[id].t_subtype_cd8_vote)))
    }) throw new IllegalStateException(name+' readback failed')
    receipt.text = gson.toJson([image:name,entry:entry.getID(),detections:actual.size(),changed:changed,classes:actual.countBy { it.getPathClass()?.toString() },polygons:geometry])
    println 'Verified '+name+' '+changed+' changes'
    System.gc()
}
def retired = ['Hematopoietic CD45⁺','Hematopoietic CD45⁻','Unclassified'] as Set
def classes = project.getPathClasses().findAll { !retired.contains(it.toString()) || usedClasses.contains(it.toString()) }
if (usedClasses.contains('T Cell')) classes += PathClass.getInstance('T Cell')
if (policy) {
    def retiredClasses = (policy.retired_detection_classes ?: []) as Set
    classes = classes.findAll { !it.toString().startsWith('Tissue cell - ') && !retiredClasses.contains(it.toString()) }
    def displays = usedClasses.findAll { it != null && allowedDisplays.contains(it) }
    displays.each { classes += catalog[it] ?: PathClass.getInstance(it) }
}
project.setPathClasses(classes.unique())
project.syncChanges()
