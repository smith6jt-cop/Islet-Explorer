/* Add cells in new Tissue minus old Tissue; never run the destructive full detection step.
 * The runner appends the usual footer, which saves only after V3_RESULT passes.
 * args[0]: augmentation config JSON. All original objects are checked before the footer saves.
 */
import com.google.gson.Gson
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.roi.GeometryTools
import qupath.lib.analysis.features.ObjectMeasurements
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.geom.prep.PreparedGeometryFactory

def cfg = new Gson().fromJson(new File(args[0]).getText('UTF-8'), Map)
def invariant = evaluate(new File(cfg.invariants_script))
def hierarchy = getCurrentHierarchy()
def oldCells = new ArrayList(hierarchy.getDetectionObjects())
def tissue = hierarchy.getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Tissue' }
if (tissue.size() != 1 || oldCells.isEmpty())
    throw new IllegalStateException('Augmentation requires one Tissue annotation and existing cells')
def imported = PathIO.readObjects(new File(cfg.outline).toPath())
if (imported.size() != 1 || imported[0].getPathClass()?.toString() != 'Tissue')
    throw new IllegalStateException('Expected exactly one Tissue in the new outline')
def target = imported[0].getROI().getGeometry()
def previous = tissue[0].getROI().getGeometry()
if (!target.isValid() || target.isEmpty()) throw new IllegalStateException('Invalid/empty new Tissue')
def prepared = PreparedGeometryFactory.prepare(target)
// Exact geometry preservation, not merely centroids. Exclusions and pLN may legitimately be outside Tissue.
def mustCover = oldCells + hierarchy.getAnnotationObjects().findAll {
    it.getPathClass()?.toString() in ['Islet', 'EndocrineCluster', 'EndocrineSingle']
}
def outside = mustCover.findAll { it.getROI() != null && it.getROI().getGeometry().getArea() > 0 &&
    !prepared.covers(it.getROI().getGeometry()) }
def outsideNuclei = oldCells.findAll { it.isCell() && it.getNucleusROI() != null &&
    it.getNucleusROI().getGeometry().getArea() > 0 &&
    !prepared.covers(it.getNucleusROI().getGeometry()) }
// A pre-existing degenerate nucleus at the origin (T1D_6534 b30903ec-...) has no area to contain.
// Preserve its original serialized ROI exactly; never create Tissue at that placeholder location.
def zeroAreaNuclei = oldCells.findAll { it.isCell() && it.getNucleusROI() != null &&
    it.getNucleusROI().getGeometry().getArea() == 0 }.collect { it.getID().toString() }
if (outside) throw new IllegalStateException("New Tissue does not cover ${outside.size()} existing geometries: " +
    outside.take(5).collect { it.getID().toString() })
if (outsideNuclei) throw new IllegalStateException("New Tissue does not cover ${outsideNuclei.size()} existing nuclei: " +
    outsideNuclei.take(5).collect { it.getID().toString() })
def saved = invariant.snapshot(hierarchy)
def inventory = new File(cfg.inventory)
if (inventory.exists() && invariant.readJSON(inventory) != saved)
    throw new IllegalStateException('Existing augmentation inventory does not match the starting image')
invariant.writeJSON(inventory, saved)
def delta = target.difference(previous)
if (!delta.isValid()) throw new IllegalStateException('New-minus-old Tissue is invalid')
def oldIDs = saved.cells.keySet() as Set
def index = new STRtree()
oldCells.each { index.insert(it.getROI().getGeometry().getEnvelopeInternal(), it) }
index.build()
def accepted = [], rejected = [outside_added: 0, overlap_old: 0, overlap_new: 0, duplicate_old: 0, duplicate_new: 0]
def result = null
def temporary = null
try {
    if (!delta.isEmpty() && delta.getArea() > 0) {
        temporary = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(delta, tissue[0].getROI().getImagePlane()),
            getPathClass('Augmentation temporary'))
        // Add directly under root; never resolve the hierarchy around the temporary parent.
        hierarchy.getRootObject().addChildObject(temporary)
        hierarchy.fireHierarchyChangedEvent(this)
        def inputChannels = ['DAPI', 'CD68', 'CD44', 'VIM', 'ECAD', 'LGALS3', 'panCK', 'CD3e', 'BActin']
        def present = getCurrentServer().getMetadata().getChannels().collect { it.getName() }
        if (!present.containsAll(inputChannels)) throw new IllegalStateException('Missing InstanSeg input channels')
        def cal = getCurrentServer().getPixelCalibration()
        if (Math.abs(cal.getPixelWidthMicrons() - 0.3774) > 1e-8 || Math.abs(cal.getPixelHeightMicrons() - 0.3774) > 1e-8)
            throw new IllegalStateException('Unexpected pixel calibration; refusing to change existing image calibration')
        hierarchy.getSelectionModel().setSelectedObjects([temporary], temporary)
        result = qupath.ext.instanseg.core.InstanSeg.builder()
            .modelPath('/home/smith6jt/QuPath/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1')
            .device('gpu0').inputChannels(inputChannels.collect { ColorTransforms.createChannelExtractor(it) })
            .tileDims(512).interTilePadding(32).nThreads(12).makeMeasurements(false)
            .randomColors(false).outputType('Default').build().detectObjects(getCurrentImageData(), [temporary])
        if (result.wasInterrupted() || result.nTilesFailed() > 0)
            throw new IllegalStateException("Incomplete InstanSeg result: ${result}")
        def candidates = new ArrayList(hierarchy.getDetectionObjects().findAll { !oldIDs.contains(it.getID().toString()) })
        // Deterministic preference for larger new candidates; existing cells always win.
        candidates.sort { a, b -> b.getROI().getArea() <=> a.getROI().getArea() ?: a.getID().toString() <=> b.getID().toString() }
        def deltaPrepared = PreparedGeometryFactory.prepare(delta)
        def factory = GeometryTools.getDefaultFactory()
        def newIndex = new org.locationtech.jts.index.quadtree.Quadtree()
        double radius = 1.0 / 0.3774
        def duplicate = { a, b ->
            def ra = a.getROI(), rb = b.getROI()
            double aa = ra.getArea(), ab = rb.getArea()
            Math.min(aa, ab) * 0.3774 * 0.3774 >= 5.0 &&
                Math.hypot(ra.getCentroidX() - rb.getCentroidX(), ra.getCentroidY() - rb.getCentroidY()) <= radius &&
                Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < 0.15
        }
        for (candidate in candidates) {
            def roi = candidate.getROI(), geom = roi.getGeometry()
            def point = factory.createPoint(new org.locationtech.jts.geom.Coordinate(roi.getCentroidX(), roi.getCentroidY()))
            if (!deltaPrepared.covers(point) || !prepared.covers(geom)) { rejected.outside_added++; continue }
            def envelope = new org.locationtech.jts.geom.Envelope(geom.getEnvelopeInternal())
            envelope.expandBy(radius)
            def neighbors = index.query(envelope)
            if (neighbors.any { geom.intersection(it.getROI().getGeometry()).getArea() > 0 }) { rejected.overlap_old++; continue }
            if (neighbors.any { duplicate(candidate, it) }) { rejected.duplicate_old++; continue }
            def nearbyNew = newIndex.query(envelope)
            if (nearbyNew.any { geom.intersection(it.getROI().getGeometry()).getArea() > 0 }) { rejected.overlap_new++; continue }
            if (nearbyNew.any { duplicate(candidate, it) }) { rejected.duplicate_new++; continue }
            accepted << candidate
            newIndex.insert(geom.getEnvelopeInternal(), candidate)
        }
        // Detach all candidates with the temporary annotation; retain only accepted objects.
        candidates.each { it.getParent()?.removeChildObject(it) }
    }
} finally {
    if (temporary != null) hierarchy.removeObject(temporary, false)
    hierarchy.getSelectionModel().clearSelection()
    hierarchy.fireHierarchyChangedEvent(this)
}
tissue[0].setROI(GeometryTools.geometryToROI(target, tissue[0].getROI().getImagePlane()))
// QuPath setROI clears this annotation's old measurements. Replace its geometry-derived shape
// fields; old Tissue-level measurements would refer to a different area and must not be retained.
ObjectMeasurements.addShapeMeasurements(tissue[0], getCurrentServer().getPixelCalibration(),
    ObjectMeasurements.ShapeFeatures.AREA, ObjectMeasurements.ShapeFeatures.LENGTH)
tissue[0].addChildObjects(accepted)
hierarchy.fireHierarchyChangedEvent(this)
def checks = invariant.verify(hierarchy, saved)
invariant.writeJSON(new File(cfg.new_ids), accepted.collect { it.getID().toString() }.sort())
invariant.requireSaveWindow(cfg)
V3_RESULT = [ok: true, changed: true, message: "Preserved ${oldCells.size()} cells; added ${accepted.size()}",
    details: [preservation: checks, added_area_um2: delta.getArea() * 0.3774 * 0.3774, new_cells: accepted.size(),
              preserved_zero_area_nuclei: zeroAreaNuclei,
              rejected: rejected, instanseg: result == null ? null : [objects: result.nObjectsDetected(),
              tiles: result.nTilesProcessed(), failed_tiles: result.nTilesFailed()],
              outputs: [cfg.inventory, cfg.new_ids]]]
