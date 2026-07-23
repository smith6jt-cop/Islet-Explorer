// export_cells_geojson.groovy — CODEX (QuPath 0.7)
//
// Export full-resolution per-cell boundaries as GeoJSON for REDSEA (phenocycler.redsea).
//
// The existing data/json/*.geojson in Islet-Explorer are ANNOTATION features (Islet, Nerve,
// Capillary, Lymphatic — 25,392 for donor 0112) and are useless here: REDSEA needs one polygon
// per DETECTION (188,371 cells for that same donor).
//
// REDSEA rasterizes these polygons into an int32 instance mask and joins the corrected per-cell
// means back on QuPath's detection UUID, so the GeoJSON MUST:
//   * be in FULL-RESOLUTION PIXEL coordinates (not microns, no downsample), and
//   * carry each feature's id == the detection's PathObject.getID() (a UUID).
// QuPath's exporter writes both by default.
//
// Output name must match redsea.resolve_paths: cells__<tag>.geojson, where <tag> is the QuPath
// image name with [^A-Za-z0-9._-] -> '_'. CODEX image names are bare zero-padded donor ids, so
// the tag is just the name: 0112 -> cells__0112.geojson.
//
// Run headlessly per image via export_cells_geojson.sh, or directly:
//   /home/smith6jt/QuPath/bin/QuPath script \
//     --project /data/QuPath_Projects/QuPath_project_IsletExplorer/project.qpproj \
//     --image 0112 --args /data/codex_phenocycler/redsea_scratch/geojson \
//     export_cells_geojson.groovy

import static qupath.lib.scripting.QP.*

// ---- config ----------------------------------------------------------------
// Output dir comes from --args (QuPath 0.7 exposes it as the `args` array); falls back to the
// default. NOTE: QuPath 0.7's `script` subcommand rejects bare -D system properties, so --args
// is the only supported channel.
// Lives on /data: /home is ~97% full and per-donor GeoJSON runs 0.9-6 GB.
def argv = this.binding.hasVariable('args') ? (this.binding.getVariable('args') as List) : []
def outPath = (argv && argv[0]) ? String.valueOf(argv[0])
                                : "/data/codex_phenocycler/redsea_scratch/geojson"
def outDir = new File(outPath)
// ----------------------------------------------------------------------------

outDir.mkdirs()
def entry = getProjectEntry()
def imageName = entry != null ? entry.getImageName()
                              : getCurrentImageData().getServer().getMetadata().getName()
def tag = imageName.replaceAll("[^A-Za-z0-9._-]", "_")
def outFile = new File(outDir, "cells__${tag}.geojson")

def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "ERROR: no detections in ${imageName} — the project entry has no cell objects"
    return
}

// Idempotent: skip a complete export so the batch driver can be re-run after a partial failure.
if (outFile.exists() && outFile.length() > 0) {
    println "SKIP ${imageName}: ${outFile} already exists (${outFile.length()} bytes)"
    return
}

// EXCLUDE_MEASUREMENTS keeps the file to geometry + id; the intensities come from the parquet,
// joined on object_id, so duplicating them here would only inflate a multi-GB file.
exportObjectsToGeoJson(detections, outFile.getAbsolutePath(),
                       "FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")

println "Wrote ${detections.size()} cell polygons -> ${outFile}"
