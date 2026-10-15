// Geometry rules adapted from IO60panc2nd/scripts/1_panc_objects.groovy.
// Pure geometry: no image metadata, source hierarchy or cell objects are changed.
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.util.AffineTransformation
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.union.UnaryUnionOp
import qupath.lib.roi.GeometryTools

def config = [method:'io60_size_adaptive_v1', edge_margin_um:5.0d,
    small_area_um2:2000.0d, medium_area_um2:10000.0d,
    small_erosion_um:1.0d, medium_erosion_um:2.0d, large_erosion_um:3.0d,
    dilation_factor:0.85d, minimum_fragment_um2:50.0d,
    hole_policy:'fill enclosed holes, as QuPath maxHoleSizeMicrons=-1; reclip to valid tissue',
    dilation_collision_policy:'Exclude shared restored area from both candidates; retain it in the raw foreground record',
    distance_units:'micrometres after independent X/Y calibration',
    minimum_eroded_fragment_source_px2:1.0d]

def polygons = { Geometry g ->
    if (g == null || g.isEmpty()) return []
    def result = []
    PolygonExtracter.getPolygons(g, result)
    result.findAll{!it.isEmpty() && it.area > 0d}.sort{a,b ->
        int c = Double.compare(a.centroid.y, b.centroid.y)
        if (!c) c = Double.compare(a.centroid.x, b.centroid.x)
        if (!c) c = Double.compare(a.area, b.area)
        if (!c) c = a.norm().toText() <=> b.norm().toText()
        c
    }
}

// A polygon intersection may also contain a detached touching line/point.
// Those zero-area pieces do not remove tissue, and legacy JTS difference
// rejects a mixed GeometryCollection. Keep only polygonal collision area.
def polygonalUnion = { List<Geometry> pieces ->
    def areaPieces=pieces.collectMany{polygons(it)}
    areaPieces.isEmpty() ? GeometryTools.createRectangle(0d,0d,0d,0d).factory.createPolygon() : UnaryUnionOp.union(areaPieces)
}

def refine = { List<Geometry> inputs, Geometry tissue, int width, int height, double sx, double sy ->
    if (![sx,sy].every{Double.isFinite(it) && it > 0d} || width <= 0 || height <= 0)
        throw new IllegalArgumentException('Positive, finite image calibration and dimensions are required')
    if (tissue == null || tissue.isEmpty() || !tissue.isValid())
        throw new IllegalArgumentException('Valid, nonempty tissue geometry is required')
    if (inputs.any{it == null || !it.isValid() || it.dimension != 2})
        throw new IllegalArgumentException('Input foreground must contain valid polygonal geometry')
    def toUm = AffineTransformation.scaleInstance(sx, sy)
    def toPixels = AffineTransformation.scaleInstance(1d/sx, 1d/sy)
    double margin = config.edge_margin_um
    if (width*sx <= 2*margin || height*sy <= 2*margin)
        throw new IllegalArgumentException('Image is smaller than the edge exclusion')
    def inner = GeometryTools.createRectangle(margin,margin,width*sx-2*margin,height*sy-2*margin)
    def allowed = toUm.transform(tissue).intersection(inner)
    def raw = inputs.collect{toUm.transform(it)}
    // Native QuPath refinement matches the source script's 50 um2 / fill-all
    // convention. Union filled shells so nested foreground islands cannot
    // become overlapping annotations. Never use a convex hull or gap closing.
    def filtered = raw.collect{GeometryTools.refineAreas(it,config.minimum_fragment_um2 as double,Double.POSITIVE_INFINITY)}
        .findAll{it != null && !it.isEmpty()}
    def envelope = filtered.isEmpty() ? allowed.factory.createPolygon() : UnaryUnionOp.union(filtered)
    def initial = polygons(envelope.intersection(allowed))
    def rawIndex = new STRtree()
    raw.eachWithIndex{g,i ->rawIndex.insert(g.envelopeInternal,i)}
    rawIndex.build()
    def output = [], audit = [], removed = []
    initial.eachWithIndex{g, parentIndex ->
        double area = g.area
        String category = area < config.small_area_um2 ? 'small' : (area < config.medium_area_um2 ? 'medium' : 'large')
        double erosion = config[category+'_erosion_um']
        double dilation = erosion * config.dilation_factor
        def eroded = g.buffer(-erosion)
        def parts = polygons(eroded).findAll{it.area >= config.minimum_eroded_fragment_source_px2*sx*sy}
        def record = [parent_index:parentIndex+1, size_category:category, input_area_um2:area,
            erosion_um:erosion, dilation_um:dilation, eroded_parts:parts.size(),
            output_parts:0, output_area_um2:0d, foreground_area_um2:0d,
            filled_background_area_um2:0d, rejected_after_dilation:0,
            status:parts.isEmpty() ? 'eroded_away' : 'retained']
        if (parts.isEmpty()) removed.add([geometry:toPixels.transform(g),reason:'eroded_away',parent_index:parentIndex+1])
        parts.each{part ->
            def expanded = part.buffer(dilation)
            def cleaned = GeometryTools.refineAreas(expanded,config.minimum_fragment_um2 as double,Double.POSITIVE_INFINITY)
            // Hole filling must never fill a tissue exclusion; dilating must
            // never undo the image-edge exclusion. Reapply both constraints.
            def clipped = cleaned.intersection(allowed)
            def candidates = polygons(clipped)
            if (candidates.isEmpty()) record.rejected_after_dilation++
            candidates.each{candidate ->
                if (candidate.area < config.minimum_fragment_um2) {
                    record.rejected_after_dilation++
                    removed.add([geometry:toPixels.transform(candidate),reason:'below_50_um2_after_refinement',parent_index:parentIndex+1])
                } else {
                    def supportPieces = rawIndex.query(candidate.envelopeInternal).collect{i ->
                        raw[i as int].intersection(candidate)
                    }.findAll{!it.isEmpty()}
                    double supportedArea = supportPieces.isEmpty() ? 0d : UnaryUnionOp.union(supportPieces).area
                    def converted = toPixels.transform(candidate)
                    if (!converted.isValid()) throw new IllegalStateException('Refinement produced invalid geometry')
                    output.add([geometry:converted,parent_index:parentIndex+1,size_category:category,
                        input_area_um2:area,erosion_um:erosion,dilation_um:dilation,
                        area_um2:candidate.area,foreground_area_um2:supportedArea,
                        filled_background_area_um2:Math.max(0d,candidate.area-supportedArea),
                        foreground_fraction:Math.min(1d,supportedArea/candidate.area)])
                    record.output_parts++
                    record.output_area_um2 += candidate.area
                    record.foreground_area_um2 += supportedArea
                    record.filled_background_area_um2 += Math.max(0d,candidate.area-supportedArea)
                }
            }
        }
        if (!parts.isEmpty() && record.output_parts == 0) record.status='below_50_um2_after_refinement'
        audit.add(record)
    }
    // Near-critical necks can split under erosion and collide again under
    // partial dilation. Remove only the shared restored area from BOTH
    // candidates, with no arbitrary assignment to the first object. This
    // preserves separation and makes the excluded area explicit.
    def proposedIndex = new STRtree()
    output.eachWithIndex{r,i ->proposedIndex.insert(r.geometry.envelopeInternal,i)}
    proposedIndex.build()
    def conflicts=[:].withDefault{[]}, shared=[]
    output.eachWithIndex{r,i ->
        proposedIndex.query(r.geometry.envelopeInternal).findAll{(it as int)<i}.each{j ->
            def overlap=r.geometry.intersection(output[j as int].geometry)
            if(overlap.area*sx*sy>1e-8) {
                conflicts[i].add(overlap);conflicts[j as int].add(overlap);shared.add(overlap)
            }
        }
    }
    def resolved=[]
    output.eachWithIndex{r,i ->
        def collision=conflicts[i]
        if(collision.isEmpty()) {
            resolved.add(r)
        } else {
            def cut=polygonalUnion(collision)
            audit[r.parent_index-1].overlap_trimmed_um2=(audit[r.parent_index-1].overlap_trimmed_um2 ?: 0d)+cut.area*sx*sy
            polygons(r.geometry.difference(cut)).each{piece ->
                if(piece.area*sx*sy < config.minimum_fragment_um2) {
                    audit[r.parent_index-1].rejected_after_dilation++
                    removed.add([geometry:piece,reason:'below_50_um2_after_overlap_exclusion',parent_index:r.parent_index])
                } else {
                    def physical=toUm.transform(piece)
                    def supported=rawIndex.query(physical.envelopeInternal).collect{j ->raw[j as int].intersection(physical)}.findAll{!it.isEmpty()}
                    double supportArea=supported.isEmpty() ? 0d : UnaryUnionOp.union(supported).area
                    def updated=new LinkedHashMap(r)
                    updated.putAll([geometry:piece,area_um2:physical.area,foreground_area_um2:supportArea,
                        filled_background_area_um2:Math.max(0d,physical.area-supportArea),
                        foreground_fraction:Math.min(1d,supportArea/physical.area)])
                    resolved.add(updated)
                }
            }
        }
    }
    output=resolved
    audit.each{r ->
        def children=output.findAll{it.parent_index==r.parent_index}
        r.output_parts=children.size()
        r.output_area_um2=children.sum{it.area_um2} ?: 0d
        r.foreground_area_um2=children.sum{it.foreground_area_um2} ?: 0d
        r.filled_background_area_um2=children.sum{it.filled_background_area_um2} ?: 0d
        r.overlap_trimmed_um2=r.overlap_trimmed_um2 ?: 0d
        if(r.output_parts==0 && r.status=='retained') r.status='below_50_um2_after_refinement'
    }
    // Reject any unresolved positive-area overlap.
    def index = new STRtree()
    output.eachWithIndex{r,i ->index.insert(r.geometry.envelopeInternal,i)}
    index.build()
    output.eachWithIndex{r,i ->
        index.query(r.geometry.envelopeInternal).findAll{(it as int) < i}.each{j ->
            double overlapArea=r.geometry.intersection(output[j as int].geometry).area*sx*sy
            if (overlapArea > 1e-6)
                throw new IllegalStateException("Refinement produced overlapping candidate regions: ${i}/${j}, parents ${r.parent_index}/${output[j as int].parent_index}, overlap ${overlapArea} um2")
        }
    }
    def summary = [input_foreground_regions:inputs.size(),
        raw_regions_removed_by_initial_fragment_filter:inputs.size()-filtered.size(),
        regions_after_filling_merging_and_clipping:initial.size(),
        parents_by_size:audit.countBy{it.size_category},
        eroded_away:audit.count{it.status=='eroded_away'},
        split_parents:audit.count{it.eroded_parts>1},
        parents_rejected_after_refinement:audit.count{it.status=='below_50_um2_after_refinement'},
        output_regions:output.size(),output_regions_by_parent_size:output.countBy{it.size_category},
        output_area_um2:output.sum{it.area_um2} ?: 0d,
        foreground_area_within_output_um2:output.sum{it.foreground_area_um2} ?: 0d,
        filled_background_area_within_output_um2:output.sum{it.filled_background_area_um2} ?: 0d,
        dilation_overlap_pairs:shared.size(),
        dilation_collision_area_excluded_um2:shared.isEmpty() ? 0d : UnaryUnionOp.union(shared).area*sx*sy,
        valid_geometry:true,positive_area_overlaps:0]
    [config:config,records:output,audit:audit,removed:removed,summary:summary]
}

return [config:config,refine:refine,polygonalUnion:polygonalUnion]
