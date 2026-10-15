// Build an independent, editable native hierarchy from schema-4 canonical inputs.
// Args: explicit donor job JSON. Original hierarchies, UUIDs, ROIs and raw means
// are preserved; corrected means and new typing are added under separate names.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.images.servers.ImageServerBuilder
import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.objects.classes.PathClass
import java.awt.image.BufferedImage
import java.nio.charset.StandardCharsets
import java.security.MessageDigest
import java.util.zip.GZIPInputStream

def argv=this.binding.getVariable('args') as List
if(argv.size()!=1)throw new IllegalArgumentException('Provide donor job JSON')
def gson=GsonTools.getInstance(true)
def jobFile=new File(argv[0].toString()).canonicalFile
def job=gson.fromJson(jobFile.text,Map.class)
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def source=new File(job.source_native.toString()).canonicalFile
def target=new File(job.output_data.toString()).canonicalFile
if(target.exists())throw new IllegalStateException('Preserve existing delivered/user-edited hierarchy')
if(sha(source)!=job.source_sha256)throw new IllegalStateException('Native source changed')
def table=new File(job.cell_table.toString())
if(sha(table)!=job.cell_table_sha256)throw new IllegalStateException('Typing table changed')
def geometry=new File(job.geometry.toString())
if(sha(geometry)!=job.geometry_sha256)throw new IllegalStateException('Geometry changed')
def original=PathIO.readHierarchy(source)
def cells=original.getCellObjects().sort{it.getID().toString()}
if(cells.size()!=job.cells.intValue())throw new IllegalStateException('Raw cell universe differs')
def originalMeasureKeys=cells[0].getMeasurements().keySet().toSet()
def digestCells={objects,keys ->
    def d=MessageDigest.getInstance('SHA-256')
    objects.sort{it.getID().toString()}.each{cell ->
        def means=new TreeMap();cell.getMeasurements().each{k,v ->if(keys.contains(k))means[k]=v}
        def r=[id:cell.getID().toString(),roi:cell.getROI().getGeometry().toText(),
               nucleus:cell.getNucleusROI()?.getGeometry()?.toText(),measurements:means]
        d.update(GsonTools.getInstance().toJson(r).getBytes(StandardCharsets.UTF_8));d.update((byte)10)
    }
    d.digest().encodeHex().toString()
}
// Every original measurement key participates, including sparse per-cell keys.
cells.each{originalMeasureKeys.addAll(it.getMeasurements().keySet())}
String before=digestCells(cells,originalMeasureKeys)
def byID=cells.collectEntries{[(it.getID().toString()):it]};def seen=new HashSet()
table.withReader('UTF-8'){reader ->
    def header=reader.readLine().split('\t',-1).toList();String line
    while((line=reader.readLine())!=null){
        def row=line.split('\t',-1);def id=row[header.indexOf('native_uuid')]
        def cell=byID[id]
        if(cell==null || !seen.add(id))throw new IllegalStateException('Missing/duplicate native cell '+id)
        String status=row[header.indexOf('classification_status')]
        String label=row[header.indexOf('base_type')]
        String lineage=row[header.indexOf('lineage')]
        String display=status=='validated' ? label : status=='lineage_supported' ? 'Lineage supported: '+lineage : status.replace('_',' ')
        cell.getMetadata().put('Previous classification',cell.getPathClass()?.toString() ?: '')
        cell.getMetadata().put('Schema4 classification status',status)
        cell.getMetadata().put('Schema4 lineage',lineage)
        cell.getMetadata().put('Schema4 phenotype state',row[header.indexOf('phenotype_state')])
        cell.setPathClass(PathClass.fromString(display))
        header.eachWithIndex{column,k ->
            if(column.startsWith('Corrected: ') || column.startsWith('Marker positive: ') || column=='Typing: endocrine member'){
                if(row[k]!='')cell.getMeasurementList().put(column,Double.parseDouble(row[k]))
            }
        }
        cell.getMeasurementList().close()
    }
}
if(seen.size()!=cells.size())throw new IllegalStateException('Some raw cells missing typing disposition')
def annotations=PathIO.readObjects(geometry)
def cores=annotations.findAll{it.getPathClass()?.toString()=='Islet'}.collectEntries{[(it.getName()):it]}
if(cores.size()!=job.islets.intValue())throw new IllegalStateException('Core identity mismatch')
def geometryText=geometry.name.endsWith('.gz') ? new GZIPInputStream(new FileInputStream(geometry)).getText('UTF-8') : geometry.text
def geometryDocument=gson.fromJson(geometryText,Map.class)
geometryDocument.features.each{feature ->
    def properties=feature.get('properties')
    def core=cores[properties.name]
    if(core!=null){
        core.getMetadata().put('Boundary assessment',properties.annotation_status?.toString() ?: 'Not recorded')
        if(properties.primary_assessable==false){
            core.setColor(255,80,80)
            core.getMetadata().put('Primary analysis','Excluded; candidate retained for inspection')
        }
    }
}
annotations.each{it.setLocked(false)}
def hierarchy=new PathObjectHierarchy()
// Clear original parentage in this in-memory copy, retaining every native cell.
hierarchy.getRootObject().addChildObjects(cells)
hierarchy.getRootObject().addChildObjects(annotations)
new File(job.core_memberships.toString()).withReader('UTF-8'){reader ->
    reader.readLine();String line
    while((line=reader.readLine())!=null){
        def fields=line.split(',',-1);String id=fields[1].split(':',2)[1]
        if(byID[id]==null || cores[fields[2]]==null)throw new IllegalStateException('Unknown exact core membership')
        cores[fields[2]].addChildObject(byID[id])
    }
}
def server=ImageServerProvider.buildServer(job.image.toString(),BufferedImage.class)
def imageData=new ImageData<BufferedImage>(server,hierarchy,ImageData.ImageType.FLUORESCENCE)
if(job.source_pixel_size_um!=null){
    double px=job.source_pixel_size_um.doubleValue()
    imageData.updateServerMetadata(new ImageServerMetadata.Builder(server.getMetadata()).pixelSizeMicrons(px,px).build())
}
imageData.setProperty('Schema4 candidate',job.candidate_root)
imageData.setProperty('Schema4 cell typing','Reviewed hierarchical all31-marker identities; raw means retained, REDSEA corrected means separately named')
imageData.setProperty('Schema4 geometry','All21 endocrine-pixel boundaries with native exact centroid memberships; flagged candidates remain reviewable')
target.parentFile.mkdirs()
PathIO.writeImageData(target,imageData)
new File(target.parentFile,'server.json').text=gson.toJson(server.getBuilder(),ImageServerBuilder.ServerBuilder.class)
def metadata=server.getMetadata();imageData.close()
def saved=PathIO.readHierarchy(target)
String after=digestCells(saved.getCellObjects(),originalMeasureKeys)
if(before!=after || saved.getCellObjects().size()!=cells.size())throw new IllegalStateException('Raw cell geometry/measurements changed on native roundtrip')
if(sha(source)!=job.source_sha256)throw new IllegalStateException('Source native data modified')
def summary=[imageType:'FLUORESCENCE',timestamp:System.currentTimeMillis(),serverMetadata:metadata,
             hierarchy:[nObjects:saved.nObjects()+1,objectTypeCounts:[Cells:cells.size(),Annotations:annotations.size()],
                        annotationClassificationCounts:saved.getAnnotationObjects().countBy{it.getPathClass()?.toString() ?: 'Unclassified'}]]
new File(target.parentFile,'summary.json').text=gson.toJson(summary)
def result=[case_id:job.case_id,status:'native_roundtrip_verified',cells:cells.size(),islets:cores.size(),
            source_data_modified:false,raw_cells_and_measurements_preserved:true,raw_digest_before:before,raw_digest_after:after,
            output_data_sha256:sha(target),job_sha256:sha(jobFile),source_native_sha256:job.source_sha256,
            geometry_sha256:job.geometry_sha256,cell_table_sha256:job.cell_table_sha256]
new File(target.parentFile,'verification.json').text=gson.toJson(result)
println "${job.case_id}: editable native schema4 hierarchy verified; ${cells.size()} cells and ${cores.size()} cores"
