// Import pixel candidates into a fresh independent review copy, preserving every cell.
// Args: fresh review project directory, prediction directory, optional case ID.
// Independent donor imports serialize shared metadata through a project lock.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.objects.classes.PathClass
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.nio.charset.StandardCharsets
import java.security.MessageDigest

def argv=this.binding.getVariable('args') as List
if(!(argv.size() in [2,3])) throw new IllegalArgumentException('Provide fresh review project and prediction directories, optionally case ID')
def root=new File(argv[0].toString()).canonicalFile
def predictionRoot=new File(argv[1].toString()).canonicalFile
def gson=GsonTools.getInstance(true)
String selectedCase=argv.size()==3 ? argv[2].toString() : null
def resultFile=selectedCase==null ? new File(root,'pixel_review_verification.json') : new File(root,'verification/'+selectedCase+'-pixel_review_verification.json')
if(resultFile.exists()) throw new IllegalStateException('Preserve delivered/user-edited review project')
def copied=gson.fromJson(new File(root,'review_manifest.json').text,Map.class)
def prediction=gson.fromJson(new File(predictionRoot,'prediction_manifest.json').text,Map.class)
if(!(prediction.status in ['refined_pixel_prediction_unreviewed','unanimous_vote_sensitivity_not_selected_primary','size_adaptive_pixel_candidates_unreviewed']))
    throw new IllegalStateException('Prediction is not a supported diagnostic review')
if(prediction.boundary_policy?.active=='io60_size_adaptive_v1' && prediction.boundary_policy.stage!='complete')
    throw new IllegalStateException('Run build_size_adaptive_islets.groovy before importing current boundary predictions')
def vectors=gson.fromJson(new File(predictionRoot,'vectorization_manifest.json').text,Map.class)
if(selectedCase==null && copied.donors.size()!=1) throw new IllegalStateException('Select a case explicitly for a multi-donor project')
def donors=copied.donors.findAll{it.case_id==(selectedCase ?: prediction.case_id)}
if(donors.size()!=1 || donors[0].case_id!=prediction.case_id) throw new IllegalStateException('Expected one matching donor')
def donor=donors[0]
def target=new File(donor.review_data.toString()).canonicalFile
def source=new File(donor.source_data.toString()).canonicalFile
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
if(!target.toPath().startsWith(root.toPath()) || Files.isSameFile(source.toPath(),target.toPath()))
    throw new IllegalStateException('Must be an independent review copy')
if(sha(source)!=donor.data_sha256 || sha(target)!=donor.data_sha256) throw new IllegalStateException('Source or review copy changed; preserve edits')
if(new File(donor.source_image.toString()).canonicalPath!=new File(prediction.source_image.toString()).canonicalPath)
    throw new IllegalStateException('Image identity differs')
def overlay=new File(predictionRoot,'review_candidates.geojson')
def record=vectors.outputs.find{it.kind=='boundary_regions'}
if(record.geometry_valid!=true || sha(overlay)!=record.sha256) throw new IllegalStateException('Native vector verification failed')
if(vectors.source_prediction_manifest_sha256!=sha(new File(predictionRoot,'prediction_manifest.json')))
    throw new IllegalStateException('Prediction metadata changed after vectorization')
def hashCells={hierarchy ->
    def d=MessageDigest.getInstance('SHA-256')
    hierarchy.getCellObjects().sort{it.getID().toString()}.each{cell ->
        def r=[id:cell.getID().toString(),parent:cell.getParent()?.getID()?.toString(),classification:cell.getPathClass()?.toString(),
            roi:cell.getROI().getGeometry().toText(),nucleus:cell.getNucleusROI()?.getGeometry()?.toText(),measurements:new TreeMap(cell.getMeasurements())]
        d.update(GsonTools.getInstance().toJson(r).getBytes(StandardCharsets.UTF_8));d.update((byte)10)
    }
    d.digest().encodeHex().toString()
}
def imageData=PathIO.readImageData(target)
def hierarchy=imageData.getHierarchy()
String before=hashCells(hierarchy)
def original=hierarchy.getAnnotationObjects().collectEntries{[(it.getID().toString()):it.getROI().getGeometry().toText()]}
hierarchy.getAnnotationObjects().each{a ->
    if(a.getPathClass()?.toString() in ['Islet','IsletExpanded']) {
        a.setPathClass(PathClass.fromString('Previous segmentation: '+a.getPathClass().toString()))
        a.setColor(150,150,150)
    }
}
def objects=PathIO.readObjects(overlay)
if(objects.size()!=record.objects.intValue() || objects.isEmpty()) throw new IllegalStateException('Candidate count mismatch')
objects.each{a ->a.setLocked(false)}
hierarchy.getRootObject().addChildObjects(objects)
hierarchy.fireHierarchyChangedEvent(this)
imageData.setProperty('Endocrine pixel review','Cyan: unreviewed endocrine pixel candidates. Grey: previous islet geometry. Existing cell parents and regional measurements refer to previous geometry. No analysis acceptance.')
def staged=new File(target.parentFile,'pixel_review.qpdata')
PathIO.writeImageData(staged,imageData);imageData.close()
imageData=null;hierarchy=null;System.gc()
def saved=PathIO.readHierarchy(staged)
String after=hashCells(saved)
if(before!=after || saved.getCellObjects().size()!=donor.cell_count.intValue()) throw new IllegalStateException('Raw cell records changed')
def annotations=saved.getAnnotationObjects().collectEntries{[(it.getID().toString()):it]}
original.each{id,wkt ->if(annotations[id]?.getROI()?.getGeometry()?.toText()!=wkt)throw new IllegalStateException('Original geometry changed')}
objects.each{a ->if(annotations[a.getID().toString()]?.getROI()?.getGeometry()?.toText()!=a.getROI().getGeometry().toText())throw new IllegalStateException('New geometry changed on serialization')}
if(sha(source)!=donor.data_sha256) throw new IllegalStateException('Processing source changed')
Files.move(staged.toPath(),target.toPath(),StandardCopyOption.REPLACE_EXISTING)
def summaryFile=new File(target.parentFile,'summary.json')
def summary=gson.fromJson(summaryFile.text,Map.class)
summary.hierarchy.annotationClassificationCounts=saved.getAnnotationObjects().countBy{it.getPathClass()?.toString() ?: 'Unclassified'}
summary.hierarchy.objectTypeCounts.Annotations=saved.getAnnotationObjects().size()
summary.hierarchy.nObjects=saved.nObjects()+1;summary.timestamp=System.currentTimeMillis()
summaryFile.text=gson.toJson(summary)
def projectFile=new File(root,'project.qpproj')
String savedProjectHash
def metadataFile=new RandomAccessFile(new File(root,'.import-metadata.lock'),'rw')
def metadataChannel=metadataFile.getChannel()
def metadataLock=metadataChannel.lock()
try {
    def project=gson.fromJson(projectFile.text,Map.class)
    def entry=project.images.find{it.entryID.intValue()==donor.entry_id.intValue()}
    if(entry==null) throw new IllegalStateException('Missing matching project image')
    entry.metadata['Boundary review']=prediction.status=='size_adaptive_pixel_candidates_unreviewed' ?
        'IO60 size-adaptive endocrine pixel candidates; original15 reference training; unreviewed' :
        'Endocrine pixel candidates; unreviewed; original15 reference training'
    entry.metadata['Delivery status']='Awaiting donor-specific visual review; not accepted replacement islets'
    projectFile.text=gson.toJson(project)
    def paletteFile=new File(root,'classifiers/classes.json')
    def palette=gson.fromJson(paletteFile.text,Map.class)
    ['Endocrine pixel candidate':-16711681,'Previous segmentation: Islet':-6908266,'Previous segmentation: IsletExpanded':-6908266].each{name,color ->
        if(!palette.pathClasses.any{it.name==name})palette.pathClasses.add([name:name,color:color])
    }
    paletteFile.text=gson.toJson(palette)
    savedProjectHash=sha(projectFile)
} finally {
    metadataLock.release();metadataChannel.close();metadataFile.close()
}
resultFile.parentFile.mkdirs()
resultFile.text=gson.toJson([status:'native_review_roundtrip_passed',case_id:prediction.case_id,
    cells:saved.getCellObjects().size(),new_boundary_candidates:objects.size(),
    cell_record_sha256_before:before,cell_record_sha256_after:after,cells_and_measurements_unchanged:true,
    source_data_unchanged:true,original_annotation_geometry_unchanged:true,new_coordinates_preserved:true,
    source_data_sha256:donor.data_sha256,review_data_sha256:sha(target),overlay_sha256:sha(overlay),
    project_sha256_at_import:savedProjectHash,model_sha256:prediction.model_sha256,accepted_islets:0,
    delivery_status:'diagnostic_review_only',operating_point:prediction.operating_point,
    prior_user_review_projects_modified:false,current_app_or_analysis_modified:false])
println "Native pixel review saved and checked: ${objects.size()} candidate regions, ${saved.getCellObjects().size()} unchanged cells"
