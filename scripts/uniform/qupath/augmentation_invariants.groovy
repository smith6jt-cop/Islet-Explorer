/* Shared read-only invariants for the additive Tissue update. Returns closures to its caller.
 * Parent links, signed distances and Tissue-only shape measurements are derived and may change;
 * all original cell raw measurements,
 * UUIDs, classes, names, cell/nucleus ROIs, and non-Tissue annotation ROIs must be preserved.
 */
import com.google.gson.Gson
import org.locationtech.jts.io.WKBWriter

def fingerprint = { obj, boolean omitROI ->
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    def put = { value ->
        byte[] bytes = String.valueOf(value).getBytes('UTF-8')
        digest.update(java.nio.ByteBuffer.allocate(4).putInt(bytes.length).array())
        digest.update(bytes)
    }
    put(obj.getName()); put(obj.getPathClass()?.toString()); put(obj.getClass().getName())
    put(obj.isLocked()); put(new TreeMap(obj.getMetadata()).toString())
    def writeROI = { roi ->
        if (roi == null) { put('NO_ROI'); return }
        put(roi.getImagePlane().toString())
        digest.update(new WKBWriter().write(roi.getGeometry()))
    }
    if (!omitROI) writeROI(obj.getROI())
    if (obj.isCell()) writeROI(obj.getNucleusROI())
    def ml = obj.getMeasurementList()
    (omitROI ? [] : ml.getNames().findAll { !it.startsWith('Signed distance') }).sort().each { name ->
        put(name); put(Long.toHexString(Double.doubleToLongBits(ml.get(name))))
    }
    digest.digest().encodeHex().toString()
}
def snapshot = { hierarchy ->
    def cells = new TreeMap(), annotations = new TreeMap()
    hierarchy.getDetectionObjects().each { cells[it.getID().toString()] = fingerprint(it, false) }
    hierarchy.getAnnotationObjects().each {
        annotations[it.getID().toString()] = fingerprint(it, it.getPathClass()?.toString() == 'Tissue')
    }
    [cells: cells, annotations: annotations]
}
def verify = { hierarchy, saved ->
    def now = snapshot(hierarchy)
    def missing = saved.cells.keySet().findAll { !now.cells.containsKey(it) }
    def changed = saved.cells.keySet().findAll { now.cells[it] != saved.cells[it] }
    if (missing || changed) throw new IllegalStateException(
        "Original cells changed: missing=${missing.size()}, changed=${changed.size()}, examples=${changed.take(5)}")
    if (now.annotations != saved.annotations)
        throw new IllegalStateException('Original annotation IDs/classes/names/measurements/ROIs changed or annotations added')
    [old_cells: saved.cells.size(), current_cells: now.cells.size(), annotations: now.annotations.size()]
}
def writeJSON = { File target, value ->
    target.getParentFile().mkdirs()
    def part = new File(target.getPath() + '.part')
    part.setText(new Gson().toJson(value), 'UTF-8')
    java.nio.file.Files.move(part.toPath(), target.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING)
}
[snapshot: snapshot, verify: verify, writeJSON: writeJSON,
 requireSaveWindow: { cfg ->
    if (!cfg.deadline_epoch_ms || System.currentTimeMillis() + 600000L >= (cfg.deadline_epoch_ms as long))
        throw new IllegalStateException('Refusing save with less than ten minutes before the driver deadline')
 },
 readJSON: { File file -> new Gson().fromJson(file.getText('UTF-8'), Object) }]
