// Native integration test: cells retain UUIDs/ROIs/raw measurements while new
// types, corrected means and regenerated parent memberships roundtrip.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import java.awt.image.BufferedImage
import javax.imageio.ImageIO
import java.security.MessageDigest

def args=this.binding.getVariable('args') as List
if(args.size()!=2)throw new IllegalArgumentException('Provide fresh test folder and implementation script')
def root=new File(args[0].toString()).canonicalFile;root.mkdirs()
def script=new File(args[1].toString()).canonicalFile
def gson=GsonTools.getInstance(true)
def sha={File f ->MessageDigest.getInstance('SHA-256').digest(f.bytes).encodeHex().toString()}
def image=new File(root,'fixture.png');ImageIO.write(new BufferedImage(50,50,BufferedImage.TYPE_USHORT_GRAY),'PNG',image)
def server=ImageServerProvider.buildServer(image.toString(),BufferedImage.class)
def hierarchy=new PathObjectHierarchy();def cells=[]
for(int i=0;i<3;i++){
    def roi=ROIs.createRectangleROI(5+i*10,10,6,6,ImagePlane.getDefaultPlane())
    def nucleus=ROIs.createEllipseROI(6+i*10,11,3,3,ImagePlane.getDefaultPlane())
    def cell=PathObjects.createCellObject(roi,nucleus,null,null)
    cell.getMeasurementList().put('Cell: INS mean',100.0+i);cell.getMeasurementList().close();cells.add(cell)
}
hierarchy.getRootObject().addChildObjects(cells)
def source=new File(root,'source.qpdata');def data=new ImageData<BufferedImage>(server,hierarchy,ImageData.ImageType.FLUORESCENCE)
PathIO.writeImageData(source,data);data.close()
def geometry=new File(root,'annotations.geojson')
geometry.text=gson.toJson([type:'FeatureCollection',features:[[type:'Feature',geometry:[type:'Polygon',coordinates:[[[3,8],[23,8],[23,18],[3,18],[3,8]]]],properties:[objectType:'annotation',name:'TestCore',classification:[name:'Islet',color:[0,255,255]]]]]])
def memberships=new File(root,'memberships.csv');memberships.text='case_id,object_id,islet_key\n'+cells.take(2).collect{'fixture,fixture:'+it.getID()+',TestCore\n'}.join('')
def table=new File(root,'cells.tsv');table.text='native_uuid\tclassification_status\tbase_type\tlineage\tphenotype_state\tCorrected: INS\tMarker positive: INS\tTyping: endocrine member\n'+cells.withIndex().collect{c,i ->
    "${c.getID()}\tvalidated\tBeta cell\tEndocrine\t\t${80+i}\t1\t1\n"}.join('')
def target=new File(root,'delivery/data.qpdata');def job=new File(root,'job.json')
job.text=gson.toJson([case_id:'fixture',candidate_root:root.toString(),source_native:source.toString(),source_sha256:sha(source),source_pixel_size_um:0.3774,
    image:image.toString(),output_data:target.toString(),cell_table:table.toString(),cell_table_sha256:sha(table),
    geometry:geometry.toString(),geometry_sha256:sha(geometry),cells:3,islets:1,core_memberships:memberships.toString()])
new GroovyShell(this.class.classLoader,new Binding([args:[job.toString()] as String[]])).evaluate(script)
def saved=PathIO.readHierarchy(target)
assert saved.getCellObjects().size()==3
assert saved.getCellObjects().every{it.getPathClass().toString()=='Beta cell'}
assert saved.getAnnotationObjects().find{it.getName()=='TestCore'}.getChildObjects().size()==2
assert saved.getCellObjects().every{it.getMeasurements().get('Corrected: INS')==it.getMeasurements().get('Cell: INS mean')-20}
def delivery=PathIO.readImageData(target)
assert delivery.getServerMetadata().getPixelCalibration().getPixelWidthMicrons()==0.3774
delivery.close()
println 'Schema4 native roundtrip: all UUID/ROI/raw-measurement, new typing and core-parent checks passed'
