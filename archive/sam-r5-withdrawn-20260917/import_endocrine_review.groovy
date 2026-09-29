// Modify only independent, unedited review copies; preserve raw cells and old geometry.
// QuPath script --args /new/review/root --args /candidate/root import_endocrine_review.groovy
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.objects.classes.PathClass
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.security.MessageDigest
import java.nio.charset.StandardCharsets

def argv = this.binding.hasVariable('args') ? this.binding.getVariable('args') as List : []
if (!(argv.size() in [2, 3, 4])) throw new IllegalArgumentException('Provide fresh review and candidate directories, optionally comma-separated donor IDs and a prior cell verification project')
def root = new File(argv[0].toString()).canonicalFile
def candidateRoot = new File(argv[1].toString()).canonicalFile
def gson = GsonTools.getInstance(true)
[root, candidateRoot].each { directory ->
    def qualityFile = new File(directory, 'quality_status.json')
    if (qualityFile.exists()) {
        def status = gson.fromJson(qualityFile.text, Map.class).status.toString()
        if (status.startsWith('withdrawn') || status.startsWith('failed'))
            throw new IllegalStateException("Boundary quality withdrawal blocks import: ${qualityFile}")
    }
}
def manifest = gson.fromJson(new File(root, 'review_manifest.json').text, Map.class)
def evidence = new File(root, 'regeneration-verification')
evidence.mkdirs()
def selection = argv.size() >= 3 && argv[2].toString() != '*' ? argv[2].toString().split(',').toSet() : null
def referenceRoot = argv.size() == 4 ? new File(argv[3].toString()).canonicalFile : null
def referenceManifest = referenceRoot == null ? null : gson.fromJson(new File(referenceRoot, 'review_manifest.json').text, Map.class)
def sha = { File file ->
    def digest = MessageDigest.getInstance('SHA-256')
    file.withInputStream { stream ->
        byte[] buffer = new byte[1024*1024]
        int n
        while ((n = stream.read(buffer)) > 0) digest.update(buffer, 0, n)
    }
    digest.digest().encodeHex().toString()
}
def hashCells = { hierarchy ->
    def digest = MessageDigest.getInstance('SHA-256')
    hierarchy.getCellObjects().sort { it.getID().toString() }.each { cell ->
        // Include UUIDs, exact cell/nuclear geometry, original parents and all measurements.
        def record = [id:cell.getID().toString(), parent:cell.getParent()?.getID()?.toString(),
            classification:cell.getPathClass()?.toString(),
            roi:cell.getROI().getGeometry().toText(),
            nucleus:cell.getNucleusROI()?.getGeometry()?.toText(),
            measurements:new TreeMap(cell.getMeasurements())]
        digest.update(GsonTools.getInstance().toJson(record).getBytes(StandardCharsets.UTF_8))
        digest.update((byte)10)
    }
    digest.digest().encodeHex().toString()
}
def reports = []
manifest.donors.each { donor ->
    if (selection != null && !(donor.case_id.toString() in selection)) return
    if (new File(evidence, "${donor.case_id}-native.json").exists())
        throw new IllegalStateException('Never rerun over a delivered donor')
    def dataFile = new File(donor.review_data.toString()).canonicalFile
    def sourceFile = new File(donor.source_data.toString()).canonicalFile
    if (!dataFile.toPath().startsWith(root.toPath()) || dataFile == sourceFile || Files.isSameFile(dataFile.toPath(), sourceFile.toPath()))
        throw new IllegalStateException('Refuse to modify source objects')
    if (sha(dataFile) != donor.data_sha256 || sha(sourceFile) != donor.data_sha256)
        throw new IllegalStateException('Review or source data have changed; preserve all edits')
    def folder = new File(candidateRoot, donor.case_id.toString())
    def summary = gson.fromJson(new File(folder, 'candidate_summary.json').text, Map.class)
    def overlayFile = new File(folder, 'review_candidates.geojson')
    if (sha(overlayFile) != summary.overlay_sha256) throw new IllegalStateException('Overlay checksum mismatch')
    println "Loading independent native copy: ${donor.case_id}"
    def imageData = PathIO.readImageData(dataFile)
    def hierarchy = imageData.getHierarchy()
    def reference = null
    String before
    if (referenceRoot != null) {
        def priorDonor = referenceManifest.donors.find { it.case_id == donor.case_id }
        reference = gson.fromJson(new File(referenceRoot, "regeneration-verification/${donor.case_id}-native.json").text, Map.class)
        if (priorDonor?.data_sha256 != donor.data_sha256 || reference.status != 'native_roundtrip_passed' ||
            reference.cells_and_measurements_unchanged != true || reference.source_data_unchanged != true)
            throw new IllegalStateException('Reference cell digest does not describe this byte-identical source')
        before = reference.cell_record_sha256
    } else before = hashCells(hierarchy)
    def originalAnnotations = hierarchy.getAnnotationObjects().collectEntries { [(it.getID().toString()):it.getROI().getGeometry().toText()] }
    int previousCores = 0, previousBands = 0
    hierarchy.getAnnotationObjects().each { annotation ->
        def cls = annotation.getPathClass()?.toString()
        if (cls in ['Islet', 'IsletExpanded']) {
            if (cls == 'Islet') previousCores++
            else previousBands++
            annotation.setPathClass(PathClass.fromString('Previous segmentation: '+cls))
            annotation.setColor(150, 150, 150)
            annotation.getMetadata().put('Geometry status', 'Previous segmentation retained for comparison; not regenerated')
        }
    }
    def objects = PathIO.readObjects(overlayFile)
    if (objects.size() != summary.new_boundary_count.intValue()+summary.separate_point_focus_annotations.intValue())
        throw new IllegalStateException('Imported overlay object count differs')
    objects.each { object ->
        object.setLocked(false)
        object.getMetadata().put('Review status', 'Unreviewed; operational hormone thresholds are not validated phenotype gates')
        object.getMetadata().put('Measurement status', 'Seed counts recalculated; peri-islet measurements and phenotype counts not recalculated')
        object.getMetadata().put('Candidate source', overlayFile.canonicalPath)
    }
    // Add as siblings without moving existing cells away from their original parents.
    hierarchy.getRootObject().addChildObjects(objects)
    hierarchy.fireHierarchyChangedEvent(this)
    imageData.setProperty('Islet candidate review', 'New green/yellow boundaries are unreviewed. Grey previous geometry, cell parents, distances and regional measurements refer to the previous segmentation. Blue annotations are seed points, not islet boundaries.')
    def staged = new File(dataFile.parentFile, 'regenerated.qpdata')
    PathIO.writeImageData(staged, imageData)
    imageData.close()
    imageData = null; hierarchy = null
    System.gc()
    println "Verifying saved geometry and every cell: ${donor.case_id}"
    def saved = PathIO.readHierarchy(staged)
    def after = hashCells(saved)
    if (before != after || saved.getCellObjects().size() != donor.cell_count.intValue())
        throw new IllegalStateException('Raw cells, parents, marker values or boundaries changed')
    def savedAnnotations = saved.getAnnotationObjects().collectEntries { [(it.getID().toString()):it] }
    originalAnnotations.each { id, geometry ->
        if (savedAnnotations[id]?.getROI()?.getGeometry()?.toText() != geometry)
            throw new IllegalStateException('An original annotation geometry changed')
    }
    objects.each { object ->
        def readBack = savedAnnotations[object.getID().toString()]
        if (readBack?.getROI()?.getGeometry()?.toText() != object.getROI().getGeometry().toText())
            throw new IllegalStateException('New boundary coordinates changed during native serialization')
    }
    Files.move(staged.toPath(), dataFile.toPath(), StandardCopyOption.REPLACE_EXISTING)
    def report = [case_id:donor.case_id, status:'native_roundtrip_passed', cells:saved.getCellObjects().size(),
        cell_record_sha256:after, cells_and_measurements_unchanged:true, original_annotation_geometry_unchanged:true,
        previous_islets:previousCores, previous_expanded_regions:previousBands,
        new_boundary_count:summary.new_boundary_count, point_focus_annotations:summary.separate_point_focus_annotations,
        new_coordinates_preserved:true, new_data_sha256:sha(dataFile), source_data_unchanged:sha(sourceFile)==donor.data_sha256,
        source_native_sha256:donor.data_sha256,
        original_cell_digest_source:referenceRoot == null ? 'computed before this save' : referenceRoot.canonicalPath,
        saved_cells_fully_hashed_in_this_run:true,
        candidate_overlay_sha256:summary.overlay_sha256,
        annotation_classes:saved.getAnnotationObjects().countBy { it.getPathClass()?.toString() ?: 'Unclassified' }]
    // Update only the independent project's thumbnail summary, not source evidence.
    def summaryFile = new File(dataFile.parentFile, 'summary.json')
    def nativeSummary = gson.fromJson(summaryFile.text, Map.class)
    nativeSummary.timestamp = System.currentTimeMillis()
    nativeSummary.hierarchy.annotationClassificationCounts = report.annotation_classes
    nativeSummary.hierarchy.objectTypeCounts.Annotations = saved.getAnnotationObjects().size()
    nativeSummary.hierarchy.nObjects = saved.nObjects()+1
    summaryFile.text = gson.toJson(nativeSummary)
    new File(evidence, "${donor.case_id}-native.json").text = gson.toJson(report)
    reports.add(report)
    saved = null; objects = null; originalAnnotations = null; savedAnnotations = null
    System.gc()
    println "Completed ${donor.case_id}: ${summary.new_boundary_count} new review boundaries; every cell preserved"
}
def allReports = manifest.donors.collect { new File(evidence, "${it.case_id}-native.json") }.findAll { it.isFile() }.collect { gson.fromJson(it.text, Map.class) }
new File(evidence, 'native-summary.json').text = gson.toJson([status:allReports.size()==manifest.donors.size() ? 'passed' : 'partial', donors:allReports,
    checked_at_utc:java.time.Instant.now().toString(), accepted_islets:0,
    candidate_root:candidateRoot.canonicalPath, original_project_modified:false])
