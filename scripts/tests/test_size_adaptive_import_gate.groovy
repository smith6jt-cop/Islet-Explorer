// Verify required-stage enforcement without opening any native object file.
// Args: importer script, fresh report folder, delivered project.qpproj.
import qupath.lib.io.GsonTools
import qupath.lib.projects.ProjectIO
import qupath.lib.roi.GeometryTools
import java.awt.image.BufferedImage

def argv=this.binding.getVariable('args') as List
if(argv.size()!=3) throw new IllegalArgumentException('Provide importer, fresh report folder, project')
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Fresh check folder required')
out.mkdirs()
def gson=GsonTools.getInstance(true)
def checks=[]
[
    [name:'new_raw_prediction',status:'refined_pixel_prediction_unreviewed',stage:'pending_native_refinement',blocked:true],
    [name:'new_cutoff_sensitivity',status:'unanimous_vote_sensitivity_not_selected_primary',stage:'pending_native_refinement',blocked:true],
    [name:'completed_adaptive',status:'size_adaptive_pixel_candidates_unreviewed',stage:'complete',blocked:false],
    [name:'legacy_reproduction',status:'refined_pixel_prediction_unreviewed',stage:null,blocked:false]
].each{test ->
    def folder=new File(out,test.name);folder.mkdirs()
    def copy=new File(folder,'review');copy.mkdirs()
    def prediction=new File(folder,'prediction');prediction.mkdirs()
    new File(copy,'review_manifest.json').text='{}'
    def manifest=[status:test.status]
    if(test.stage!=null) manifest.boundary_policy=[active:'io60_size_adaptive_v1',stage:test.stage]
    new File(prediction,'prediction_manifest.json').text=gson.toJson(manifest)
    boolean stageBlocked=false, advancedToVectorCheck=false
    try {
        def binding=new Binding([args:[copy.canonicalPath,prediction.canonicalPath]])
        new GroovyShell(getClass().getClassLoader(),binding).evaluate(new File(argv[0].toString()))
    } catch(IllegalStateException e) {
        if(!e.message.contains('Run build_size_adaptive_islets.groovy')) throw e
        stageBlocked=true
    } catch(FileNotFoundException e) {
        if(!e.message.contains('vectorization_manifest.json')) throw e
        advancedToVectorCheck=true
    }
    if(stageBlocked!=test.blocked || advancedToVectorCheck==test.blocked)
        throw new IllegalStateException('Unexpected import stage enforcement: '+test.name)
    checks.add([test:test.name,passed:true])
}
def projectFile=new File(argv[2].toString()).canonicalFile
def project=ProjectIO.loadProject(projectFile,BufferedImage.class)
def stagingFile=new File(projectFile.parentFile,'staging_manifest.json')
def expectedCount=stagingFile.exists() ? gson.fromJson(stagingFile.text,Map.class).donors.size() : 1
if(project.getImageList().size()!=expectedCount || project.getImageList().count{it.getImageName().contains('HPAP')}!=1)
    throw new IllegalStateException('Wrong review image')
if(!project.getPathClasses().any{it.toString()=='Endocrine pixel candidate'})
    throw new IllegalStateException('Missing review class')
def uris=project.getImageList().collectMany{it.getURIs()}
if(uris.toSet().size()!=expectedCount) throw new IllegalStateException('Duplicate source images')
uris.each{uri ->
    if(uri.scheme!='file' || !new File(uri).isFile()) throw new IllegalStateException('Missing linked image')
}
new File(out,'verification.json').text=gson.toJson([status:'passed',import_gate_checks:checks,
    native_project_load:'passed',project:projectFile.canonicalPath,
    native_geometry_grid_units:1d/GeometryTools.getDefaultFactory().getPrecisionModel().getScale(),
    native_coordinate_sequence_factory:GeometryTools.getDefaultFactory().getCoordinateSequenceFactory().toString(),
    images:project.getImageList().collect{it.getImageName()},image_count:expectedCount,source_projects_modified:false])
println 'PASS four import-stage checks and native project load'
