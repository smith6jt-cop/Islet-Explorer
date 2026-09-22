// Read the original 15 donors without saving the source project or image data.
// Args: original project.qpproj, fresh reference directory.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import java.security.MessageDigest

def argv = this.binding.getVariable('args') as List
if (argv.size() != 2) throw new IllegalArgumentException('Provide source project and fresh output directory')
def projectFile = new File(argv[0].toString()).canonicalFile
def out = new File(argv[1].toString()).canonicalFile
if (out.exists()) throw new IllegalStateException('Preserve existing references')
out.mkdirs()
def gson = GsonTools.getInstance(true)
def project = gson.fromJson(projectFile.text, Map.class)
def expected = ['112','6356','6450','6479','6505','6516','6521','6533','6534','6538','6548','6549','6550','6551','6563'].toSet()
def selected = project.images.findAll { it.imageName.toString().replaceFirst('^0+', '') in expected }
if (selected.size() != 15) throw new IllegalStateException('Original 15 donor mapping is incomplete')
def sha = { File file ->
    def d = MessageDigest.getInstance('SHA-256')
    file.withInputStream { s -> byte[] b = new byte[1048576]; int n; while ((n=s.read(b))>0) d.update(b,0,n) }
    d.digest().encodeHex().toString()
}
def records = []
selected.each { entry ->
    String cid = entry.imageName.toString().replaceFirst('^0+', '')
    def nativeFile = new File(projectFile.parentFile, "data/${entry.entryID.intValue()}/data.qpdata")
    String before = sha(nativeFile)
    long modified = nativeFile.lastModified(), bytes = nativeFile.length()
    println "Reading original donor ${cid} annotations"
    def hierarchy = PathIO.readHierarchy(nativeFile)
    def all = hierarchy.getAnnotationObjects()
    def references = all.findAll { it.getPathClass()?.toString() in ['Islet','Tissue'] }
    def islets = references.findAll { it.getPathClass()?.toString() == 'Islet' }
    if (islets.isEmpty() || references.size() == islets.size()) throw new IllegalStateException('Missing reference classes')
    def donorDir = new File(out,cid); donorDir.mkdirs()
    def target = new File(donorDir,'reference_annotations.geojson')
    PathIO.exportObjectsAsGeoJSON(target,references,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
    def b = entry.serverBuilder
    def metadata = b.metadata
    while (b.containsKey('builder')) b = b.builder
    def source = new File(new URI(b.uri.toString())).canonicalFile
    if (!source.isFile()) throw new IllegalStateException("Missing original image ${source}")
    if (nativeFile.lastModified()!=modified || nativeFile.length()!=bytes) throw new IllegalStateException('Source changed during export')
    def record = [case_id:cid,entry_id:entry.entryID.intValue(),source_image:source.toString(),
        source_image_bytes:source.length(),source_image_mtime_ms:source.lastModified(),
        source_native:nativeFile.canonicalPath,source_native_sha256:before,
        source_native_bytes:bytes,source_native_mtime_ms:modified,
        source_channel_names:metadata.channels.collect{it.name},
        image_shape_cyx:[metadata.channels.size(),metadata.height.intValue(),metadata.width.intValue()],
        pixel_size_um:0.3774,reference_annotations:target.canonicalPath,reference_sha256:sha(target),
        islet_annotations:islets.size(),islets_with_SAM_quality:islets.count{it.getMeasurementList().containsKey('SAM Quality')},
        annotation_classes:all.countBy{it.getPathClass()?.toString() ?: 'Unclassified'},
        manual_review_provenance:'Not encoded in the source project; reference regions are not pixel-level ground truth',
        source_modified:false]
    new File(donorDir,'reference_manifest.json').text=gson.toJson(record)
    records.add(record)
    println "Exported ${cid}: ${islets.size()} reference islets"
    hierarchy=null;all=null;references=null;islets=null
    System.gc()
}
new File(out,'reference_manifest.json').text=gson.toJson([
    status:'source_references_exported',source_project:projectFile.canonicalPath,source_project_sha256:sha(projectFile),
    qupath_version:qupath.lib.common.GeneralTools.getVersion(),donors:records,
    target_donors_excluded_from_training:['78','6090','6232','6520','6539','6558'],
    manual_pixel_training_annotations_found:false,source_project_modified:false])
