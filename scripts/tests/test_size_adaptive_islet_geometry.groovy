// Run with QuPath's installed JVM/classes. Tests the same helper used in delivery.
import org.locationtech.jts.geom.util.AffineTransformation
import org.locationtech.jts.operation.union.UnaryUnionOp
import qupath.lib.roi.GeometryTools

try {
def argv = this.binding.getVariable('args') as List
def rules = new GroovyShell(getClass().getClassLoader()).evaluate(new File(argv[0]))
def box = {x,y,w,h ->GeometryTools.createRectangle(x as double,y as double,w as double,h as double)}
def tissue = box(0,0,1000,1000)
def run = {geometries ->rules.refine(geometries,tissue,1000,1000,1d,1d)}
def checked=[]

def sizes=run([box(30,30,30,30),box(100,100,60,60),box(250,250,110,110),
    box(450,450,40,50),box(600,600,100,100)])
assert sizes.audit*.size_category==['small','medium','large','medium','large']
assert sizes.audit*.erosion_um==[1d,2d,3d,2d,3d]
assert sizes.audit*.dilation_um==[.85d,1.7d,2.55d,1.7d,2.55d]
checked << 'small/medium/large rules and exact 2000/10000 um2 thresholds'

def dumbbell=UnaryUnionOp.union([box(100,100,30,30),box(145,100,30,30),box(129,114.5,17,1)])
def split=run([dumbbell])
assert split.records.size()==2 && split.summary.split_parents==1
assert split.records[0].geometry.distance(split.records[1].geometry)>0d
checked << 'touching islets split and remain separate after 85 percent restoration'

def nearCritical=UnaryUnionOp.union([box(100,100,30,30),box(131,100,30,30),box(129,114.005,3,1.99)])
def collision=run([nearCritical])
assert collision.records.size()==2 && collision.summary.dilation_overlap_pairs>0
assert collision.summary.dilation_collision_area_excluded_um2>0d
assert collision.records[0].geometry.intersection(collision.records[1].geometry).area < 1e-6
checked << 'near-critical necks cannot be rejoined by colliding dilation fronts'

def areaCollision=box(200,200,10,10)
def touchingLine=box(240,240,10,10).boundary
def mixedCollision=areaCollision.factory.createGeometryCollection([areaCollision,touchingLine] as org.locationtech.jts.geom.Geometry[])
def polygonalCut=rules.polygonalUnion([mixedCollision])
assert polygonalCut.area==100d && polygonalCut.dimension==2
assert box(190,190,40,40).difference(polygonalCut).area==1500d
checked << 'mixed area-and-touching-line collisions subtract only positive area'

def ring=box(100,100,40,40).difference(box(110,110,20,20))
def island=box(117,117,6,6)
def hole=run([ring,island])
assert hole.records.size()==1 && hole.records[0].filled_background_area_um2>350d
assert hole.records[0].geometry.covers(box(115,115,1,1))
checked << 'native fill-all-hole policy records unsupported fill and merges nested islands'

def cutTissue=tissue.difference(box(110,110,20,20))
def cut=rules.refine([ring,island],cutTissue,1000,1000,1d,1d)
assert cut.records.every{it.geometry.intersection(box(110,110,20,20)).area < 1e-6}
checked << 'tissue exclusions survive hole filling and dilation'

def edge=run([box(0,0,35,35)])
assert edge.records.every{it.geometry.envelopeInternal.minX>=5 && it.geometry.envelopeInternal.minY>=5}
def edgeStrip=run([box(0,20,4,100)])
assert edgeStrip.records.isEmpty()
checked << 'five micrometre image-edge clipping before and after refinement'

def lost=run([box(100,100,100,1),box(250,100,6,6),box(300,100,7.2,7.2)])
assert lost.records.isEmpty()
assert lost.summary.eroded_away==1 && lost.summary.raw_regions_removed_by_initial_fragment_filter==1
assert lost.summary.parents_rejected_after_refinement==1
checked << 'thin erased fragments and pre/post 50 um2 exclusions are accounted for'

def reference=run([dumbbell,ring])
double sx=.3774d,sy=.62d
def scale=AffineTransformation.scaleInstance(1d/sx,1d/sy)
def calibrated=rules.refine([scale.transform(dumbbell),scale.transform(ring)],scale.transform(tissue),
    Math.ceil(1000/sx) as int,Math.ceil(1000/sy) as int,sx,sy)
assert calibrated.records.size()==reference.records.size()
calibrated.records.eachWithIndex{r,i ->
    def restored=AffineTransformation.scaleInstance(sx,sy).transform(r.geometry)
    assert restored.symDifference(reference.records[i].geometry).area < .01
}
checked << 'source coordinate projection and anisotropic physical distances'

def again=run([ring,dumbbell])
assert again.records*.geometry*.toText()==reference.records*.geometry*.toText()
checked << 'determinism independent of input object order'
assert run([]).records.isEmpty()
boolean rejected=false
try {rules.refine([ring],tissue,1000,1000,Double.NaN,1d)} catch(IllegalArgumentException e){rejected=true}
assert rejected
checked << 'empty foreground and invalid calibration'
println 'PASS '+checked.size()+' native morphology checks: '+checked.join('; ')
} catch(AssertionError e) {
    // QuPath's CLI closes its worker pools on Exceptions, but not Errors.
    throw new IllegalStateException(e.message,e)
}
