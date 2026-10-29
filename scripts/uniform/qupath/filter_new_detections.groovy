/* Explicit decisions only. Runs against an independent staging project, never the raw source.
 * args[0]: config JSON. modes apply (remove/save if nonempty) and verify (fresh reload, read-only).
 */
import com.google.gson.Gson

def cfg = new Gson().fromJson(new File(args[0]).getText('UTF-8'), Map)
def helper = evaluate(new File(cfg.helper))
def hierarchy = getCurrentHierarchy()
def rejected = new HashSet(cfg.excluded_ids)
def original = new HashSet(cfg.original_ids)
def rawNew = new HashSet(cfg.raw_new_ids)
def beforeFile = new File(cfg.before_inventory)
def details
if (cfg.mode == 'apply') {
    def before = helper.remove(hierarchy, original, rawNew, rejected)
    helper.writeJSON(beforeFile, before)
    details = helper.verify(hierarchy, before, rejected)
    if (!rejected.isEmpty()) {
        hierarchy.fireHierarchyChangedEvent(this)
        // Dynamic annotation counts must be regenerated; stored annotation measurements stay unchanged.
        saveAnnotationMeasurements(cfg.annotation_output.toString())
        getProjectEntry().saveImageData(getCurrentImageData())
    }
    details.saved = !rejected.isEmpty()
} else if (cfg.mode == 'verify') {
    def before = helper.readJSON(beforeFile)
    details = helper.verify(hierarchy, before, rejected)
    details.saved = false
} else throw new IllegalArgumentException('Expected apply or verify mode')
details.mode = cfg.mode
details.image = getProjectEntry().getImageName()
details.ok = true
helper.writeJSON(new File(cfg.report), details)
println('[new_detection_filter] ' + details)
