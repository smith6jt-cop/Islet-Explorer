/*
 * import_review.groovy (QuPath 0.7): prepare one image of the review copy (data/reviews/qupath-review-20261027/)
 * so every object the current root (v15) rests on can be checked against the image. Owner request, 2026-10-27 ("I need to
 * verify all annotations, detections, and cell typing in QuPath").
 *
 * The copy already holds, unchanged from the run data/development/uniform-20donor-20261027/: the Tissue outline, the
 * InstanSeg cell detections, and the Islet / EndocrineCluster / EndocrineSingle objects. This adds:
 *   - each detection's v15 cell type (data/development/uniform-20donor-20261027-v15/cell_types/) as its PathClass,
 *     matched by UUID; cell_class, refined_class, immune_definition, cell_region and islet_name go into the
 *     detection's metadata;
 *   - each detection the analysis left out (<image>.excluded.csv, from the run's typing/calls and typing/redsea) a class
 *     "Excluded - <reason>" and its recorded reason in the metadata (excluded_reason); a detection in neither file
 *     stops the script;
 *   - the Capillary, Nerve and Lymphatic objects the run used (structures/<Class>/), locked, added below the root
 *     so no cell changes parent.
 * Rerunning replaces the structures it added. Run on the review copy only:
 *   JAVA_TOOL_OPTIONS=-Xmx96g /home/smith6jt/QuPath/bin/QuPath script -p <review>/project.qpproj -i <image> \
 *       -a <assignments directory> import_review.groovy
 */
import qupath.lib.io.PathIO
import qupath.lib.objects.classes.PathClass
import java.awt.Color

String folder = (args != null && args.size() > 0 && args[0]?.trim()) ? args[0].trim() :
    "/data/panc_CODEX/Islet-Explorer/data/reviews/qupath-review-20261027-assignments"
def entry = getProjectEntry()
String imageName = entry != null ? entry.getImageName() : getCurrentImageData().getServer().getMetadata().getName()
def csv = new File(folder, imageName + ".csv")
if (!csv.isFile()) throw new IllegalStateException("no assignments for '${imageName}' at ${csv}")

// Colors: immune types strong, endocrine types pale, exocrine and untyped faint; structures distinct.
def colors = [
    "T Cell CD8a\u207a": "#E41A1C", "T Cell CD4\u207a": "#FF7F00", "T Cell": "#F781BF", "B Cell": "#FFFF33",
    "Immune": "#1B9E77", "Endocrine Beta": "#9ECAE1", "Endocrine Alpha": "#A1D99B", "Endocrine Delta": "#BCBDDC",
    "Endocrine CHGA\u207a hormone-negative": "#8D6E63", "Endocrine": "#BDBDBD", "Exocrine": "#E0E0E0",
    "Excluded - cell QC": "#555555", "Excluded - exocrine inside islet border": "#A0522D",
    "Excluded - hormone under 2x runner-up": "#9932CC", "Excluded - single without usable hormone": "#6A5ACD",
    "Excluded - ColIV top marker": "#8B7D6B", "Excluded - lineage ambiguous": "#B8860B",
    "Excluded - doublet suspected": "#FF1493",
    "CD3e\u207a": "#C9A3A3", "CD20\u207a": "#C9C2A3", "CD68\u207a": "#A3C9C2", "CD163\u207a": "#A3B5C9",
    "ECAD\u207a": "#B5A3C9", "panCK\u207a": "#C9A3BD",
    // cells no rule claims, named by their top marker (owner, 2026-10-27): muted, one color each
    "B3TUBB\u207a": "#B8A47E", "CD31\u207a": "#C98B8B", "CD34\u207a": "#D4A5A5", "CD4\u207a": "#9FB7C9",
    "CD8a\u207a": "#8FA8D0", "GCG\u207a": "#B5C99A", "INS\u207a": "#A9C4D9", "LGALS3\u207a": "#C3B1D6",
    "PDPN\u207a": "#D9B38C", "SMA\u207a": "#A3A3C2", "SST\u207a": "#C2C2A3", "VIM\u207a": "#B7C9B7",
    "CD45\u207a (top marker)": "#8CC2B5",
    // the CD45 rule's class (was "Immune"; owner, 2026-10-27)
    "CD45\u207a": "#1B9E77",
    "Capillary": "#FF2D2D", "Nerve": "#FFD700", "Lymphatic": "#E6994D", "Collagen IV": "#D82DBD"]
// Class names are the display names (owner, 2026-10-27: broad name first, superscript plus and minus), e.g.
// "Macrophage CD68\u207aCD163\u207bLGALS3\u207aCD4\u207b"; the data label is kept in the metadata (phenotype_label).
// Macrophages: hue by CD68 / CD163 (both blue, CD68 only cyan, CD163 only violet), lighter for each of LGALS3 and CD4
// that is negative.
String PLUS = "\u207a", MINUS = "\u207b"
def macrophageColor = { String name ->
    float hue = name.contains("CD68" + PLUS + "CD163" + PLUS) ? 0.58f : (name.contains("CD68" + PLUS) ? 0.50f : 0.75f)
    int negatives = (name.contains("LGALS3" + MINUS) ? 1 : 0) + (name.contains("CD4" + MINUS) ? 1 : 0)
    Color.getHSBColor(hue, (float) (1.0f - 0.3f * negatives), (float) (1.0f - 0.15f * negatives))
}
// colors written by write_assignments.py into the project's classes.json win (final typing run, 2026-10-28)
def projectColors = [:]
def cj = getProject() != null ? new File(getProject().getPath().getParent().toFile(), "classifiers/classes.json") : null
if (cj != null && cj.isFile()) {
    com.google.gson.JsonParser.parseString(cj.getText("UTF-8")).getAsJsonObject().getAsJsonArray("pathClasses").each { e ->
        def o = e.getAsJsonObject(); projectColors[o.get("name").getAsString()] = o.get("color").getAsInt() }
}
def classOf = { String name ->
    if (projectColors.containsKey(name) && !name.startsWith("Macrophage ")) return PathClass.fromString(name, projectColors[name])
    def c = name.startsWith("Macrophage ") ? macrophageColor(name) : Color.decode(colors[name] ?: "#808080")
    PathClass.fromString(name, c.getRGB())
}

def hierarchy = getCurrentHierarchy()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }

// Cell types.
def headers = [:]
def readRows = { File file ->
    def out = [:]
    def header = null
    file.eachLine("UTF-8") { line ->
        def f = line.split(",", -1)
        if (header == null) { header = f.toList(); headers[file.getName()] = header; return }
        def rec = [:]; header.eachWithIndex { h, i -> rec[h] = i < f.size() ? f[i] : "" }
        out[rec["object_id"]] = rec
    }
    out
}
def rows = readRows(csv)
// every column but image, object_id and phenotype goes into the detection's metadata; earlier keys are cleared
def metaKeys = headers[csv.getName()] - ["image", "object_id", "phenotype"]
def oldKeys = (["phenotype_label", "cell_class", "refined_class", "immune_definition", "cell_region", "islet_name",
                "typing_route", "start_class", "compartment", "compartment_status", "lineage", "lineage_status",
                "subtype_status", "status", "named_by", "tags", "confirmed", "confirming_markers"] + metaKeys).unique()
def excludedCsv = new File(folder, imageName + ".excluded.csv")
if (!excludedCsv.isFile()) throw new IllegalStateException("missing ${excludedCsv}")
def excluded = readRows(excludedCsv)
def detections = hierarchy.getDetectionObjects()
int matched = 0, untyped = 0, left = 0
def counts = new TreeMap<String, Integer>()
detections.each { d ->
    def rec = rows[d.getID().toString()]
    String cls
    if (rec == null) {
        def ex = excluded[d.getID().toString()]
        oldKeys.each { d.getMetadata().remove(it) }
        if (ex == null) { cls = "Untyped"; untyped++ }
        else { cls = ex["class"]; left++; d.getMetadata().put("excluded_reason", ex["reason"]) }
    } else {
        d.getMetadata().remove("excluded_reason")
        cls = rec["phenotype"] ?: "Untyped"; matched++
        oldKeys.each { d.getMetadata().remove(it) }
        metaKeys.each { k ->
            def v = rec[k]; if (v != null && !v.isEmpty()) d.getMetadata().put(k, v) }
    }
    d.setPathClass(classOf(cls))
    counts[cls] = (counts[cls] ?: 0) + 1
}
println "[review] ${imageName}: detections ${detections.size()}, typed ${matched} of ${rows.size()} rows, excluded ${left} of ${excluded.size()} rows, in neither ${untyped}"
counts.each { k, v -> println "[review]   ${k}: ${v}" }
if (matched != rows.size()) throw new IllegalStateException("${rows.size() - matched} typed cells have no detection with that UUID")
if (left != excluded.size() || untyped != 0) throw new IllegalStateException("excluded rows unmatched or detections in neither file")

// Structures: replace any earlier import, then add below the root without resolving the hierarchy.
def structureClasses = ["Capillary", "Nerve", "Lymphatic", "Collagen IV"]
def parentsBefore = detections.collectEntries { [(it.getID()): it.getParent()?.getID()] }
def earlier = hierarchy.getAnnotationObjects().findAll { clsOf(it) in structureClasses }
if (!earlier.isEmpty()) hierarchy.removeObjects(earlier, true)
structureClasses.each { c ->
    def f = new File(folder, "structures/${c}__${imageName}.geojson")
    if (!f.isFile()) throw new IllegalStateException("missing ${f}")
    def objs = PathIO.readObjects(f.toPath())
    objs.each { o -> o.setPathClass(classOf(c)); o.setLocked(true); hierarchy.addObjectBelowParent(hierarchy.getRootObject(), o, false) }
    println "[review]   structures ${c}: ${objs.size()}"
}
hierarchy.fireHierarchyChangedEvent(this)
int moved = detections.count { parentsBefore[it.getID()] != it.getParent()?.getID() }
if (moved != 0) throw new IllegalStateException("${moved} detections changed parent")
def annCounts = new TreeMap<String, Integer>()
hierarchy.getAnnotationObjects().each { annCounts[clsOf(it)] = (annCounts[clsOf(it)] ?: 0) + 1 }
println "[review]   annotations: ${annCounts}"
if (entry != null) { entry.saveImageData(getCurrentImageData()); println "[review] saved ${imageName}" }
else println "[review] no project entry: not saved"
