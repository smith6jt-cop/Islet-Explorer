// Export approved PolygonROI vertices without the JTS 0.01 px rounding step.
// Args: promoted native data, fresh export folder. Never writes native data.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.roi.PolygonROI
import qupath.lib.roi.GeometryROI
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.MultiPolygon
import java.security.MessageDigest
import java.util.zip.GZIPOutputStream

def argv=this.binding.getVariable('args') as List
if(argv.size()!=2) throw new IllegalArgumentException('Provide promoted native data and fresh export folder')
def source=new File(argv[0].toString()).canonicalFile
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve previous spatial audit')
def gson=GsonTools.getInstance(true)
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
String before=sha(source)
def hierarchy=PathIO.readHierarchy(source)
def active=hierarchy.getAnnotationObjects().findAll{it.getPathClass()?.toString() in ['Islet','IsletExpanded','Tissue','Exclude','Nerve','Capillary','Lymphatic']}
def polygonCoords={g ->
    def rings=[g.getExteriorRing().getCoordinates().collect{[it.x,it.y]}]
    for(int i=0;i<g.getNumInteriorRing();i++) rings.add(g.getInteriorRingN(i).getCoordinates().collect{[it.x,it.y]})
    rings
}
def features=active.collect{object ->
    def feature=GsonTools.getInstance().fromJson(GsonTools.getInstance().toJson(object),Map.class)
    if(object.getPathClass()?.toString()=='Islet') {
        def roi=object.getROI()
        if(roi instanceof PolygonROI) {
            def coords=roi.getAllPoints().collect{[it.getX(),it.getY()]}
            coords.add(coords[0])
            feature.geometry=[type:'Polygon',coordinates:[coords]]
        } else if(roi instanceof GeometryROI) {
            def g=roi.getGeometry()
            if(g instanceof Polygon) feature.geometry=[type:'Polygon',coordinates:polygonCoords(g)]
            else if(g instanceof MultiPolygon) feature.geometry=[type:'MultiPolygon',coordinates:(0..<g.getNumGeometries()).collect{polygonCoords(g.getGeometryN(it))}]
            else throw new IllegalStateException('Unexpected area geometry '+g.getGeometryType())
        } else throw new IllegalStateException('Explicit containment/export handling required for '+roi.getClass().getName())
        feature.get('properties').put('native_roi_type',roi.getClass().getSimpleName())
    }
    feature
}
out.mkdirs()
def path=new File(out,'annotations.geojson.gz')
path.withOutputStream{s ->new GZIPOutputStream(s).withWriter('UTF-8'){w ->
    GsonTools.getInstance().toJson([type:'FeatureCollection',features:features],w)
}}
def cell=hierarchy.getCellObjects().find{it.getID().toString()=='f24597cd-146a-4cd8-8be5-8426218ca3d2'}
def core=active.find{it.getName()=='AdaptiveIslet_78_00247' && it.getPathClass()?.toString()=='Islet'}
def diagnostic=[object_uuid:cell.getID().toString(),islet_key:core.getName(),
    x_source_px:cell.getROI().getCentroidX(),y_source_px:cell.getROI().getCentroidY(),
    native_contains:core.getROI().contains(cell.getROI().getCentroidX(),cell.getROI().getCentroidY()),
    roi_type:core.getROI().getClass().getSimpleName()]
if(sha(source)!=before) throw new IllegalStateException('Native source changed')
def record=[status:'exact_native_polygon_vertices_exported',native_data:source.canonicalPath,native_data_sha256:before,
    geometry_sha256:sha(path),cores:active.count{it.getPathClass()?.toString()=='Islet'},
    core_roi_types:active.findAll{it.getPathClass()?.toString()=='Islet'}.countBy{it.getROI().getClass().getSimpleName()},
    core_roi_type_by_name:active.findAll{it.getPathClass()?.toString()=='Islet'}.collectEntries{[(it.getName()):it.getROI().getClass().getSimpleName()]},
    core_coordinates:'PolygonROI getAllPoints and GeometryROI exact stored coordinates; no export-grid rounding',
    other_coordinates:'Preserved QuPath native GeoJSON export',native_data_unmodified:true,boundary_diagnostic:diagnostic]
new File(out,'verification.json').text=gson.toJson(record)
println gson.toJson(record)
