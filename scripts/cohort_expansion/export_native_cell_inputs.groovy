// Read-only native hierarchy export, without opening or loading the image pixels.
// Args: data.qpdata, fresh output directory. UUIDs are the canonical join key.
import qupath.lib.io.PathIO
import qupath.lib.io.GsonTools
import java.security.MessageDigest

def argv=this.binding.getVariable('args') as List
if(!(argv.size() in [2,3])) throw new IllegalArgumentException('Provide native data, fresh output directory, optional centroids-only')
def centroidsOnly=argv.size()==3 && argv[2].toString()=='centroids-only'
def source=new File(argv[0].toString()).canonicalFile
def out=new File(argv[1].toString()).canonicalFile
if(out.exists()) throw new IllegalStateException('Output exists; preserve prior exports')
def sha={File f ->
    def d=MessageDigest.getInstance('SHA-256')
    f.withInputStream{s ->byte[] b=new byte[1048576];int n;while((n=s.read(b))>0)d.update(b,0,n)}
    d.digest().encodeHex().toString()
}
def before=sha(source)
def hierarchy=PathIO.readHierarchy(source)
def cells=hierarchy.getCellObjects()
if(cells.isEmpty()) throw new IllegalStateException('Native hierarchy contains no cells')
if(cells.collect{it.getID().toString()}.toSet().size()!=cells.size()) throw new IllegalStateException('Duplicate cell UUID')
out.mkdirs()
def geometry=new File(out,'cells.geojson')
if(!centroidsOnly) PathIO.exportObjectsAsGeoJSON(geometry,cells,PathIO.GeoJsonExportOptions.FEATURE_COLLECTION,PathIO.GeoJsonExportOptions.EXCLUDE_MEASUREMENTS)
def centroids=new File(out,'native_centroids.tsv')
centroids.withWriter('UTF-8'){w ->
    w.write('object_id\tx_source_px\ty_source_px\tcell_area_px2\tnucleus_area_px2\n')
    cells.each{cell ->
        def r=cell.getROI();def n=cell.getNucleusROI()
        w.write([cell.getID().toString(),r.getCentroidX(),r.getCentroidY(),r.getArea(),n==null ? '' : n.getArea()].join('\t')+'\n')
    }
}
if(sha(source)!=before) throw new IllegalStateException('Native input changed during export')
new File(out,'manifest.json').text=GsonTools.getInstance(true).toJson([
    status:'native_cell_inputs_exported',source_native:source.canonicalPath,source_sha256:before,
    cells:cells.size(),geometry_sha256:centroidsOnly ? null : sha(geometry),centroid_sha256:sha(centroids),
    coordinates:'full resolution source pixels',source_modified:false,
    qupath_version:qupath.lib.common.GeneralTools.getVersion()])
println "Exported ${cells.size()} native cells from ${source}"
