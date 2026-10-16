/*
 * CODEX PORT (Islet-Explorer, 2026-10-15): copied from Islet-Explorer-Senior's IO60 v3 pipeline
 * (scripts/uniform/PORTED_FROM.json records the source file and checksum). Only the CODEX
 * parameters below were changed, and on 2026-10-16 the owner changed which cells an outline takes (below).
 * WHISKER_UM keeps Senior's 0.63475 um: the strips it removes are steps of InstanSeg's outlines, which InstanSeg draws
 * on its model's own 0.5 um grid in both pipelines (CODEX images are downsampled 0.5 / 0.3774 for it), so they are
 * 0.5 or 1.0 um wide here as in Senior. Everything else is the Senior script.
 */
/*
 * 1b_islet_outlines_to_cells.groovy -- runs after script 1 (islet outlines) on an image that has its cells (script 0).
 *
 * The maintainer, 2026-10-02, after looking at the per-image classifier's outlines:
 *
 *   "The islet outlines overall look very good, but make the following changes: All outlines need to erode while
 *    stopping at nearest cell segmentation outline. Distinction between hormone single-cell and new class "Hormone
 *    Fragment" by cell segmentation outline: if one whole cell fits in islet outline, use cell segmentation; if not,
 *    label as hormone fragment, the object erodes to marker outline and markers do not contribute to cell type."
 *
 * The owner, 2026-10-16, after looking at the result: "The endocrine objects outlines need to also pick up the cell
 * segmentations. There should not be multiple annotations in one cell outline; among other problems this causes a
 * single-cell to be counted as a cluster. There should never be an annotation for endocrine object that is smaller than
 * a cell." That replaces the Hormone Fragment of 2026-10-02 and Senior's fit tolerance.
 *
 * What it does to every Islet annotation script 1 made (geometry only; no measurement of any cell is read):
 *
 *   1. The outline PICKS UP every detection it covers by at least PICK_SHARE of the detection's area; a detection two
 *      outlines cover goes to the one covering more of it. An outline that picks up no cell that way (hormone signal
 *      inside one cell, or spread across a few) picks up the cell holding most of it, unless another outline already
 *      has that cell; largest overlaps first, so two outlines inside one cell give it one object. Every other detection
 *      the outline touches is CUT BY it. A cell = a detection with a nucleus.
 *   2. The outline takes the cells it picks up whole, and loses the ones it cuts: the detections it cuts are taken
 *      out of it, so its edge runs along their outlines; the strips left between two neighbouring cut cells are
 *      removed (WHISKER_UM: InstanSeg's outlines step in its own pixels, so those strips are 0.985 or 1.97 image pixels
 *      wide; everything under 2.5 px goes); the picked-up detections are added whole, also where they reach beyond
 *      the outline.
 *   3. What is left is split into its connected pieces; pieces holding parts of one cell (InstanSeg draws some cells in
 *      parts) stay one piece. A piece holding two or more cells is an object with that
 *      piece as its outline. A piece holding exactly one cell is an object whose outline IS that cell's segmentation
 *      outline. A piece holding none is dropped.
 *   4. An outline that picked up no cell (hormone signal outside every cell) is removed. No object is smaller than a
 *      cell, and no cell lies in two objects (check_islet_objects.py checks both on the exports).
 *
 *   5. Every object left is then DESIGNATED (the maintainer, 2026-10-02: "why are there no designations of
 *      EndocrineCluster or EndocrineSingle?" -- script 1 knows the one class Islet, and the sizes went with the
 *      detector that was retired on 10-01):
 *        EndocrineSingle    exactly one whole cell (their "hormone single-cell"; its outline is that cell's outline)
 *        EndocrineCluster   two or more whole cells, under ISLET_MIN_DIAMETER_UM across
 *        Islet              two or more whole cells, ISLET_MIN_DIAMETER_UM across or more (Drotar 2025: 35 um)
 *      "Across" is the mean of QuPath's own Max diameter and Min diameter of the outline (minimum bounding circle,
 *      minimum caliper width) -- Drotar's mean of the largest and smallest diameter; it is written on the object as
 *      'Mean caliper diameter um'. A single is decided by the cell, not by size: a third of the one-cell objects
 *      are over 10 um across, and a size rule called such a cell a cluster. Labels only: no outline changes.
 *
 * Every object this step has finished carries the measurement 'Whole cells' (how many cells lie whole in it). That is
 * also how a re-run recognises its own output: it erodes only Islets without it, and designates the rest again (so a
 * re-run after ISLET_MIN_DIAMETER_UM changed re-labels without touching an outline). What an earlier run made from
 * outlines script 1 has since drawn again -- EndocrineCluster, EndocrineSingle and Hormone Fragment objects; script 1
 * removes its own class Islet only -- is removed first: that is the case when no finished Islet is left and most of
 * those objects lie under a new outline.
 *
 * islet_outlines_to_cells.py is the same rule in Python; islet_classifier_check.py compares the two.
 */

import org.locationtech.jts.algorithm.MinimumBoundingCircle
import org.locationtech.jts.algorithm.MinimumDiameter
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.PrecisionModel
import org.locationtech.jts.geom.TopologyException
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import org.locationtech.jts.geom.util.GeometryFixer
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.buffer.BufferOp
import org.locationtech.jts.operation.buffer.BufferParameters
import org.locationtech.jts.operation.overlayng.OverlayNG
import org.locationtech.jts.operation.overlayng.OverlayNGRobust
import org.locationtech.jts.operation.overlayng.UnaryUnionNG
import qupath.lib.images.ImageData
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.roi.GeometryTools
import qupath.opencv.ml.pixel.PixelClassifierTools

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
double PICK_SHARE = 0.5               // owner 2026-10-16: an outline covering at least this share of a cell picks the whole cell up
double WHISKER_UM = 0.63475           // strips of an eroded outline narrower than twice this (2.5 px) are removed
String CLASS_FRAGMENT = 'Hormone Fragment'   // no longer made (owner 2026-10-16); an earlier run's are removed
String CLASS_CLUSTER = 'EndocrineCluster'
String CLASS_SINGLE = 'EndocrineSingle'
double ISLET_MIN_DIAMETER_UM = 35.0   // two or more whole cells and this far across or more: an Islet (Drotar 2025)
String KEY_CELLS = 'Whole cells'
// ====================================================

def imageData = getCurrentImageData()
if (imageData == null) throw new IllegalStateException('[01b] No image open.')
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
double pxUm = Math.sqrt(cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons())
double whiskerPx = WHISKER_UM / pxUm
String tag = getProjectEntry().getImageName().replaceAll('[^A-Za-z0-9._-]+', '_')
long t0 = System.currentTimeMillis()

def detections = new ArrayList(getDetectionObjects())
if (detections.isEmpty()) throw new IllegalStateException('[01b] No detections: the cells (script 0) must exist before islet outlines can end on them.')
def clsOf = { o -> o.getPathClass()?.toString() }
def ENDOCRINE = [CLASS_ISLET, CLASS_CLUSTER, CLASS_SINGLE]
String KEY_DIAMETER = "Mean caliper diameter ${cal.getPixelWidthUnit()}".toString()
def isDone = { o -> clsOf(o) in ENDOCRINE && o.getROI() != null && o.getROI().isArea() && !Double.isNaN(getMeas(o, KEY_CELLS)) }
def fresh = getAnnotationObjects().findAll { clsOf(it) == CLASS_ISLET && it.getROI() != null && it.getROI().isArea() && Double.isNaN(getMeas(it, KEY_CELLS)) }
def finished = getAnnotationObjects().findAll { isDone(it) }
def oldFragments = getAnnotationObjects().findAll { clsOf(it) == CLASS_FRAGMENT }
println "[01b] ${detections.size()} detections; ${fresh.size()} Islet annotations from script 1; ${finished.size()} objects already on cell outlines " +
        "${new TreeMap(finished.countBy { clsOf(it) })}; ${oldFragments.size()} ${CLASS_FRAGMENT} annotations"

def details = [tag: tag, pick_share: PICK_SHARE, whisker_um: WHISKER_UM, islet_min_diameter_um: ISLET_MIN_DIAMETER_UM,
               detections: detections.size(), islets_from_script_1: fresh.size(), islets_already_done: finished.size()]

// ---- 5. designation: one whole cell = EndocrineSingle; more = EndocrineCluster, or Islet from ISLET_MIN_DIAMETER_UM across ----
def designate = {
    def counts = new TreeMap([(CLASS_ISLET): 0, (CLASS_CLUSTER): 0, (CLASS_SINGLE): 0])
    int relabelled = 0, measured = 0
    getAnnotationObjects().findAll { isDone(it) }.each { o ->
        def g = o.getROI().getGeometry()
        double across = (2.0d * new MinimumBoundingCircle(g).getRadius() + new MinimumDiameter(g).getLength()) / 2.0d * pxUm
        String want = getMeas(o, KEY_CELLS) < 2 ? CLASS_SINGLE : (across >= ISLET_MIN_DIAMETER_UM ? CLASS_ISLET : CLASS_CLUSTER)
        if (clsOf(o) != want) { o.setPathClass(getPathClass(want)); relabelled++ }
        double had = getMeas(o, KEY_DIAMETER)
        if (Double.isNaN(had) || Math.abs(had - across) > 1e-9) { putMeas(o, KEY_DIAMETER, across); measured++ }
        counts[want] = counts[want] + 1
    }
    return [counts: counts, relabelled: relabelled, measured: measured]
}

if (fresh.isEmpty()) {
    // nothing script 1 has made since the last run of this step: no outline is touched; the designation is made again
    def des = designate()
    if (des.relabelled > 0) fireHierarchyUpdate()
    details += [designation: des.counts, relabelled: des.relabelled]
    V3_RESULT = [ok: !finished.isEmpty() || !oldFragments.isEmpty(), changed: des.relabelled > 0 || des.measured > 0,
                 message: (finished.isEmpty() && oldFragments.isEmpty()) ? 'no Islet annotation on this image (run script 1 first)'
                         : ("no new outline from script 1: all ${finished.size()} objects already end on cell outlines; designation ${des.counts} " +
                            "(${des.relabelled} re-labelled)").toString(),
                 details: details]
    println "[01b] ${V3_RESULT.message}"
} else {
    // What an earlier run made from outlines that script 1 has since drawn again: script 1 removes its own class
    // (Islet) only, so the EndocrineCluster / EndocrineSingle / Hormone Fragment objects of that run are still here.
    // That is the case when no finished Islet is left AND most of those objects lie under a new outline. (An image
    // that never had an Islet-sized object and gets one hand-drawn outline keeps what it has: nothing lies under it.)
    int staleFragments = 0, staleDesignated = 0
    def derived = finished.findAll { clsOf(it) != CLASS_ISLET } + oldFragments
    if (!finished.any { clsOf(it) == CLASS_ISLET } && !derived.isEmpty()) {
        def freshTree = new STRtree()
        def freshPrep = new IdentityHashMap()
        fresh.each { f -> freshTree.insert(f.getROI().getGeometry().getEnvelopeInternal(), f) }
        freshTree.build()
        int under = derived.count { d ->
            def pt = d.getROI().getGeometry().getInteriorPoint()
            freshTree.query(pt.getEnvelopeInternal()).any { f ->
                def pg = freshPrep.get(f)
                if (pg == null) { pg = PreparedGeometryFactory.prepare(f.getROI().getGeometry()); freshPrep.put(f, pg) }
                pg.contains(pt)
            }
        }
        if (2 * under >= derived.size()) {
            staleFragments = oldFragments.size()
            staleDesignated = derived.size() - staleFragments
            removeObjects(derived, true)
            finished = finished.findAll { clsOf(it) == CLASS_ISLET }
            println "[01b] Removed ${staleDesignated} ${CLASS_CLUSTER} / ${CLASS_SINGLE} and ${staleFragments} ${CLASS_FRAGMENT} annotations left from an " +
                    "earlier run (${under} of ${derived.size()} lie under an outline script 1 has drawn since)"
        } else {
            println "[01b] Kept ${derived.size()} ${CLASS_CLUSTER} / ${CLASS_SINGLE} / ${CLASS_FRAGMENT} annotations of an earlier run: only ${under} lie under " +
                    "one of the ${fresh.size()} new outlines, so script 1 has not drawn their outlines again"
        }
    }

    // Overlays of pixel-staircase outlines now and then meet a geometry that is not valid (a ring touching itself);
    // JTS then throws "side location conflict". Repair the inputs and try again; last, snap to 0.001 px.
    def snap = new PrecisionModel(1000.0d)
    int nRepaired = 0
    def fixed = { Geometry g -> (g == null || g.isEmpty() || g.isValid()) ? g : GeometryFixer.fix(g) }
    // Everything here is area. An intersection of two areas that also share a stretch of boundary comes back as
    // polygons PLUS lines, and JTS refuses such a mixed collection as the input of the next overlay.
    def areaOnly = { Geometry g ->
        if (g == null || g.isEmpty() || g instanceof org.locationtech.jts.geom.Polygonal) return g
        def polys = []
        PolygonExtracter.getPolygons(g, polys)
        return polys.isEmpty() ? g.getFactory().createPolygon() : g.getFactory().buildGeometry(polys)
    }
    def overlay = { Geometry a0, Geometry b0, int op ->
        def a = areaOnly(a0), b = areaOnly(b0)
        try {
            return areaOnly(OverlayNGRobust.overlay(a, b, op))
        } catch (TopologyException e) {
            nRepaired++
            try {
                return areaOnly(OverlayNGRobust.overlay(fixed(a), fixed(b), op))
            } catch (TopologyException e2) {
                return areaOnly(OverlayNG.overlay(fixed(a), fixed(b), op, snap))
            }
        }
    }
    def unionOf = { Collection<Geometry> gs ->
        if (gs.size() == 1) return gs[0]
        try {
            return areaOnly(OverlayNGRobust.union(gs))
        } catch (TopologyException e) {
            nRepaired++
            return areaOnly(UnaryUnionNG.union(gs.collect { fixed(it) }, snap))
        }
    }

    // ---- spatial index of the detections (envelopes only; a cell's geometry is built when an outline comes near it) ----
    int n = detections.size()
    def tree = new STRtree()
    for (int i = 0; i < n; i++) {
        def r = detections[i].getROI()
        tree.insert(new Envelope(r.getBoundsX(), r.getBoundsX() + r.getBoundsWidth(), r.getBoundsY(), r.getBoundsY() + r.getBoundsHeight()), i)
    }
    tree.build()
    def geomCache = new HashMap<Integer, Geometry>()
    def cellGeom = { int i ->
        def g = geomCache.get(i)
        if (g == null) {
            g = fixed(detections[i].getROI().getGeometry())
            geomCache.put(i, g)
        }
        return g
    }
    def hasNucleus = { int i -> def d = detections[i]; d.isCell() && d.getNucleusROI() != null }
    def prepFactory = new PreparedGeometryFactory()
    // The opening that removes the strips: shrink by WHISKER, grow back. ROUND joins: with mitre joins JTS 1.20 returns
    // almost nothing for the negative buffer of some multi-part outlines (found 2026-10-02 on one outline in ~2,000:
    // 21,717 px^2 in 40 parts shrank to 4 px^2, and the islet came out as 16 scraps), while every part alone, and the
    // whole with round joins, shrink correctly. So the result is also checked: the whole cells are much wider than a
    // strip, and an opening that takes away more than a tenth of them has failed. Then: part by part; then not at all.
    def roundJoin = new BufferParameters()
    int nOpenedByPart = 0, nOpeningSkipped = 0
    def openOnce = { Geometry g ->
        def shrunk = BufferOp.bufferOp(g, -whiskerPx, roundJoin)
        if (shrunk.isEmpty()) return shrunk
        // an opening never adds anything; the arcs of the round joins are polygons, though, and can reach a few
        // thousandths of a pixel beyond the outline, so the result is cut back to it
        return overlay(areaOnly(fixed(BufferOp.bufferOp(shrunk, whiskerPx, roundJoin))), g, OverlayNG.INTERSECTION)
    }
    def keepsCells = { Geometry result, Geometry before, Geometry cellsUnion ->
        double want = overlay(before, cellsUnion, OverlayNG.INTERSECTION).getArea()
        return want <= 0 || overlay(result, cellsUnion, OverlayNG.INTERSECTION).getArea() >= 0.9 * want
    }
    def opened = { Geometry g, Geometry cellsUnion ->
        if (whiskerPx <= 0 || g.isEmpty()) return g
        def result = openOnce(g)
        if (keepsCells(result, g, cellsUnion)) return result
        def parts = []
        PolygonExtracter.getPolygons(g, parts)
        def byPart = parts.collect { openOnce(it) }.findAll { !it.isEmpty() }
        result = byPart.isEmpty() ? g.getFactory().createPolygon() : unionOf(byPart)
        if (keepsCells(result, g, cellsUnion)) { nOpenedByPart++; return result }
        nOpeningSkipped++
        return g
    }

    // ---- 1. which cells each outline picks up (owner, 2026-10-16: "the endocrine object outlines need to also pick up
    // the cell segmentations"; no cell holds two objects; no object is smaller than a cell) ----
    // a. An outline covering at least PICK_SHARE of a cell picks the whole cell up; a cell such outlines share goes to the
    //    one covering most of it. (Senior took a cell only when it also reached no farther than one classifier pixel
    //    beyond the outline, and cut the others out.)
    // b. An outline that picks up no cell with a nucleus that way -- hormone signal inside one cell, or across a few --
    //    picks up the cell with a nucleus holding most of it, if no other outline has that cell. Largest overlaps first,
    //    so two outlines inside one cell give that cell one object. (Senior made such an outline a Hormone Fragment
    //    shaped like the hormone's own pixels, smaller than a cell.)
    // c. An outline that ends up with no cell (hormone signal outside every cell) is removed.
    int nIslets = fresh.size()
    def geoms = fresh.collect { fixed(it.getROI().getGeometry()) }
    def touching = new ArrayList<List<Integer>>(nIslets)
    def overlapOf = new ArrayList<Map<Integer, Double>>(nIslets)
    def home = new HashMap<Integer, Integer>()
    def homeArea = new HashMap<Integer, Double>()
    for (int o = 0; o < nIslets; o++) {
        def g = geoms[o]
        def prepG = prepFactory.create(g)
        def mine = new ArrayList<Integer>()
        def areas = new HashMap<Integer, Double>()
        for (cand in tree.query(g.getEnvelopeInternal())) {
            int i = cand as int
            def c = cellGeom(i)
            if (!prepG.intersects(c)) continue
            double inside = overlay(g, c, OverlayNG.INTERSECTION).getArea()
            if (inside <= 0) continue                                  // touching along an edge is not overlapping
            mine << i
            areas.put(i, inside)
            if (inside >= PICK_SHARE * c.getArea()) {
                def best = homeArea.get(i)
                if (best == null || inside > best) { home.put(i, o); homeArea.put(i, inside) }
            }
        }
        touching << mine
        overlapOf << areas
    }
    def picksANucleatedCell = new boolean[nIslets]
    home.each { i, o -> if (hasNucleus(i)) picksANucleatedCell[o] = true }
    def minorityPicks = []
    for (int o = 0; o < nIslets; o++) {
        if (picksANucleatedCell[o]) continue
        touching[o].findAll { hasNucleus(it) }.each { i -> minorityPicks << [o: o, i: i, area: overlapOf[o].get(i)] }
    }
    int nPickedByMinority = 0
    minorityPicks.sort { -it.area }.each { p ->
        if (picksANucleatedCell[p.o] || home.containsKey(p.i)) return
        home.put(p.i, p.o); homeArea.put(p.i, p.area)
        picksANucleatedCell[p.o] = true
        nPickedByMinority++
    }

    // ---- 2 + 3. each outline: take out the cells it does not pick up, add the whole cells it does; split; keep the
    // pieces that hold a cell ----
    def toAdd = []
    def noCell = []
    int nSplit = 0, nSingle = 0, nMulti = 0, nPiecesDropped = 0
    long wholeCells = 0, wholeDetections = 0
    for (int o = 0; o < nIslets; o++) {
        def obj = fresh[o]
        def whole = touching[o].findAll { home.get(it) == o }
        if (!whole.any { hasNucleus(it) }) { noCell << obj; continue }
        def cut = touching[o].findAll { home.get(it) != o }
        Geometry region = geoms[o]
        def wholeUnion = unionOf(whole.collect { cellGeom(it) })
        if (!cut.isEmpty()) region = overlay(region, unionOf(cut.collect { cellGeom(it) }), OverlayNG.DIFFERENCE)
        region = opened(region, wholeUnion)
        region = overlay(region, wholeUnion, OverlayNG.UNION)
        def pieces = []
        PolygonExtracter.getPolygons(region, pieces)
        // About 7 % of InstanSeg cell outlines come in parts (6,609 of 93,727 on ND_6232, median 14 % of the cell
        // outside its largest part). Pieces holding parts of one picked-up cell are one object, so a cell is never split
        // between objects or left partly outside its object.
        int np = pieces.size()
        int[] group = new int[np]
        for (int p = 0; p < np; p++) group[p] = p
        def root = { int a -> while (group[a] != a) { group[a] = group[group[a]]; a = group[a] }; a }
        def prepPieces = pieces.collect { prepFactory.create(it) }
        def pieceOf = new HashMap<Integer, Integer>()
        for (i in whole) {
            def c = cellGeom(i)
            def hit = (0..<np).findAll { prepPieces[it].intersects(c) && overlay(pieces[it], c, OverlayNG.INTERSECTION).getArea() > 0 }
            if (hit.isEmpty()) continue
            for (int k = 1; k < hit.size(); k++) { int a = root(hit[0]), b = root(hit[k]); if (a != b) group[a] = b }
            pieceOf.put(i, hit[0])
        }
        def groups = new LinkedHashMap<Integer, List<Integer>>()
        for (int p = 0; p < np; p++) groups.computeIfAbsent(root(p), { new ArrayList<Integer>() }).add(p)
        def plane = obj.getROI().getImagePlane()
        def parent = obj.getParent()
        int kept = 0
        for (g in groups.entrySet()) {
            def piece = g.value.size() == 1 ? pieces[g.value[0]] : unionOf(g.value.collect { pieces[it] })
            def members = whole.findAll { pieceOf.containsKey(it) && root(pieceOf.get(it)) == g.key }
            def nucleated = members.findAll { hasNucleus(it) }
            if (nucleated.isEmpty()) { nPiecesDropped++; continue }
            // one whole cell: "use cell segmentation" -- the outline IS that cell's outline
            Geometry outline = nucleated.size() == 1 ? cellGeom(nucleated[0]) : piece
            if (nucleated.size() == 1) nSingle++ else nMulti++
            wholeCells += nucleated.size()
            wholeDetections += (nucleated.size() == 1 ? 1 : members.size())
            def target = obj
            if (kept > 0) {
                target = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(outline, plane), obj.getPathClass())
                toAdd << [parent: parent, object: target]
            } else {
                obj.setROI(GeometryTools.geometryToROI(outline, plane))
            }
            putMeas(target, KEY_CELLS, nucleated.size() as double)
            kept++
        }
        if (kept == 0) noCell << obj
        else if (kept > 1) nSplit++
    }
    toAdd.groupBy { it.parent }.each { parent, items ->
        if (parent != null) parent.addChildObjects(items.collect { it.object }) else addObjects(items.collect { it.object })
    }

    // ---- 4. an outline that picked up no cell is removed: no endocrine object is smaller than a cell ----
    if (!noCell.isEmpty()) removeObjects(noCell, true)
    def des = designate()
    fireHierarchyUpdate()

    int nIsletsOut = nSingle + nMulti
    details += [islets: nIsletsOut, islets_one_whole_cell: nSingle, islets_two_or_more: nMulti, outlines_split: nSplit,
                pieces_without_a_whole_cell_dropped: nPiecesDropped, whole_cells: wholeCells, whole_detections: wholeDetections,
                cells_picked_up_by_a_minority_overlap: nPickedByMinority, outlines_without_a_cell_removed: noCell.size(),
                stale_fragments_removed: staleFragments, stale_designated_removed: staleDesignated,
                designation: des.counts, relabelled: des.relabelled, overlays_repaired: nRepaired,
                openings_done_part_by_part: nOpenedByPart, openings_skipped: nOpeningSkipped,
                seconds: (System.currentTimeMillis() - t0) / 1000.0]
    boolean ok = nIsletsOut > 0
    V3_RESULT = [ok: ok, changed: true,
                 message: ("${fresh.size()} outlines from script 1 -> ${nIsletsOut} objects on cell outlines (${nSingle} one cell, ${nMulti} more; " +
                         "${nSplit} outlines split; ${nPickedByMinority} cells picked up by an outline covering less than half of them); " +
                         "${noCell.size()} outlines without a cell removed; " +
                         "designation ${des.counts}").toString(),
                 details: details]
    println "[01b] ${V3_RESULT.message} (${details.seconds} s)"
}
