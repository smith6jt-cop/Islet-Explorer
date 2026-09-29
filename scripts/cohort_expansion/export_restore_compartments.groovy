// Read-only UUID-selected native means for the locked RESTORE Membrane+Cell fit.
// Args: data.qpdata, requested UUID text file, markers JSON list, fresh output dir.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import java.security.MessageDigest

def argv=this.binding.getVariable('args') as List
if(argv.size()!=4) throw new IllegalArgumentException('Provide native data, UUID list, marker JSON, output directory')
def source=new File(argv[0].toString()).canonicalFile
def requested=new File(argv[1].toString()).readLines('UTF-8').findAll{!it.isBlank()}.toSet()
def markers=GsonTools.getInstance().fromJson(new File(argv[2].toString()).text,List)
def out=new File(argv[3].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Preserve the previous native export')
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s -> byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def before=sha(source)
def cells=PathIO.readHierarchy(source).getCellObjects().findAll{requested.contains(it.getID().toString())}
if(cells.size()!=requested.size()) throw new IllegalStateException('Requested UUIDs absent or duplicated')
def columns=['object_id']+['Membrane','Cell'].collectMany{c -> markers.collect{m -> "${c}__${m}"}}
def measurementNames=cells.collectMany{it.getMeasurementList().getNames()}.toSet()
def sourceNames=[:]
['Membrane','Cell'].each{c -> markers.each{m ->
    def variants=["${c}: ${m}: Mean".toString(),"${c}: ${m.toString().replace('PGP9.5','PGP9_5')}: Mean".toString()].unique()
    def matches=variants.findAll{measurementNames.contains(it)}
    if(matches.size()!=1) throw new IllegalStateException("Missing/ambiguous native measurement ${c} ${m}: ${matches}")
    sourceNames["${c}__${m}".toString()]=matches[0]
}}
out.mkdirs()
def table=new File(out,'compartments.tsv')
table.withWriter('UTF-8'){w ->
    w.write(columns.join('\t')+'\n')
    cells.sort{a,b -> a.getID().toString()<=>b.getID().toString()}.each{cell ->
        def ml=cell.getMeasurementList()
        def values=[cell.getID().toString()]+columns.drop(1).collect{Double.toString(ml.get(sourceNames[it]))}
        w.write(values.join('\t')+'\n')
    }
}
if(sha(source)!=before) throw new IllegalStateException('Native input changed during export')
new File(out,'manifest.json').text=GsonTools.getInstance(true).toJson([
    source_native:source.canonicalPath,source_sha256:before,cells:cells.size(),
    requested_uuid_sha256:sha(new File(argv[1].toString())),markers_sha256:sha(new File(argv[2].toString())),
    table_sha256:sha(table),measurement_names:sourceNames,
    numeric_encoding:'java_double_to_string',source_modified:false,
    qupath_version:qupath.lib.common.GeneralTools.getVersion()])
println "Exported ${cells.size()} UUID-selected native Membrane+Cell profiles"
