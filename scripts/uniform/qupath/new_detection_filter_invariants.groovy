/* Full-object invariants for deletion of explicitly rejected NEW leaf detections only.
 * This deliberately includes signed distances, parent UUIDs and Tissue ROI: none may change.
 */
import com.google.gson.Gson
import org.locationtech.jts.io.WKBWriter

def fingerprint = { obj ->
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    def put = { value ->
        byte[] bytes = String.valueOf(value).getBytes('UTF-8')
        digest.update(java.nio.ByteBuffer.allocate(4).putInt(bytes.length).array()); digest.update(bytes)
    }
    put(obj.getID()); put(obj.getClass().getName()); put(obj.getName())
    put(obj.getPathClass()?.toString()); put(obj.isLocked()); put(new TreeMap(obj.getMetadata()))
    put(obj.getColor()); put(Long.toHexString(Double.doubleToRawLongBits(obj.getClassProbability())))
    if (obj.isAnnotation()) put(obj.getDescription())
    def parent = obj.getParent()
    put(parent == null ? 'NONE' : parent.isRootObject() ? 'ROOT' : parent.getID().toString())
    def roiHash = { roi ->
        if (roi == null) { put('NO_ROI'); return }
        put(roi.getImagePlane()); digest.update(new WKBWriter().write(roi.getGeometry()))
    }
    roiHash(obj.getROI()); if (obj.isCell()) roiHash(obj.getNucleusROI())
    def ml = obj.getMeasurementList()
    new ArrayList(ml.getNames()).sort().each { name -> put(name); put(Long.toHexString(Double.doubleToRawLongBits(ml.get(name)))) }
    digest.digest().encodeHex().toString()
}
def snapshot = { hierarchy ->
    def cells = new TreeMap(), annotations = new TreeMap(), counts = new TreeMap()
    hierarchy.getDetectionObjects().each { cells[it.getID().toString()] = fingerprint(it) }
    hierarchy.getAnnotationObjects().each { annotations[it.getID().toString()] = fingerprint(it); counts[it.getID().toString()] = 0 }
    hierarchy.getDetectionObjects().each { cell ->
        def parent = cell.getParent()
        while (parent != null) {
            if (parent.isAnnotation()) counts[parent.getID().toString()] += 1
            parent = parent.getParent()
        }
    }
    [cells: cells, annotations: annotations, annotation_detection_counts: counts]
}
def verify = { hierarchy, before, rejected ->
    def now = snapshot(hierarchy)
    def expected = new TreeMap(before.cells); rejected.each { expected.remove(it) }
    if (now.cells != expected) throw new IllegalStateException('Retained native cells/parents/measurements/ROIs changed or inventory differs')
    if (now.annotations != before.annotations) throw new IllegalStateException('Native annotation metadata/measurements/parents/ROIs changed')
    [retained_cells: now.cells.size(), excluded_cells: rejected.size(), annotations: now.annotations.size(),
     retained_full_signatures_match: true, annotation_full_signatures_match: true,
     annotation_detection_counts: now.annotation_detection_counts]
}
def remove = { hierarchy, originals, rawNew, rejected ->
    def current = hierarchy.getDetectionObjects()
    def ids = current.collect { it.getID().toString() } as Set
    def expected = new HashSet(originals); expected.addAll(rawNew)
    if (!Collections.disjoint(originals, rawNew) || !rawNew.containsAll(rejected))
        throw new IllegalStateException('Rejection is not a subset of disjoint new IDs')
    if (ids.size() != expected.size() || !ids.containsAll(expected))
        throw new IllegalStateException('Native source inventory differs from original plus raw new IDs')
    def doomed = current.findAll { rejected.contains(it.getID().toString()) }
    if (doomed.any { !it.getChildObjects().isEmpty() })
        throw new IllegalStateException('Refusing to remove a detection with children')
    def before = snapshot(hierarchy)
    if (!doomed.isEmpty()) hierarchy.removeObjects(doomed, true)
    verify(hierarchy, before, rejected)
    before
}
def writeJSON = { File target, value ->
    target.getParentFile().mkdirs()
    def part = new File(target.getPath() + '.part')
    part.setText(new Gson().toJson(value), 'UTF-8')
    java.nio.file.Files.move(part.toPath(), target.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
}
[snapshot: snapshot, verify: verify, remove: remove, writeJSON: writeJSON,
 readJSON: { File file -> new Gson().fromJson(file.getText('UTF-8'), Map) }]
