"""Scientific identity and measurement contracts for the expansion importer."""
import importlib.util
import json
import gzip
from pathlib import Path
import tempfile
import unittest
import sys

import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parents[1] / "cohort_expansion/prepare_analysis.py"
spec = importlib.util.spec_from_file_location("cohort_analysis", MODULE)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)
audit_spec = importlib.util.spec_from_file_location("cohort_audit", MODULE.with_name("audit_preparation.py"))
audit = importlib.util.module_from_spec(audit_spec)
previous_analysis = sys.modules.get("prepare_analysis")
sys.modules["prepare_analysis"] = analysis
try:
    audit_spec.loader.exec_module(audit)
finally:
    if previous_analysis is None:
        sys.modules.pop("prepare_analysis", None)
    else:
        sys.modules["prepare_analysis"] = previous_analysis


class ImportTests(unittest.TestCase):
    def setUp(self):
        self.row = {"case_id": "78", "source_filename": "hp.tif", "inclusion_status": "selected"}
        self.meta = {"Age": 18., "Gender": "F", "Donor Status": "T1D", "donor_source": "HPAP"}
        self.frame = pd.DataFrame({"Image": ["hp.tif", "hp.tif"], "Object ID": ["a", "b"],
                                   "Parent": ["Islet_1", "Islet_1_exp20um"],
                                   "Centroid X µm": [3.774, 7.548], "Centroid Y µm": [7.548, 3.774],
                                   "Cell: INS: Mean": [100., 0.], "Cell: GCG: Mean": [2., 20.],
                                   "Cell: Area µm^2": [50., 60.],
                                   "Signed distance to annotation Islet µm": [-3., 5.]})

    def test_raw_import_preserves_marker_and_micrometre_coordinates(self):
        data = analysis.convert_frame(self.frame, ["GCG", "INS"], self.row, self.meta)
        np.testing.assert_array_equal(data.X, [[2., 100.], [20., 0.]])
        np.testing.assert_allclose(data.obs["X_centroid"], [3.774, 7.548])
        np.testing.assert_allclose(data.layers["log"], np.log1p(data.X))
        self.assertEqual(list(data.obs.index), ["78:a", "78:b"])
        self.assertEqual(list(data.obs["Dist to Closest Islet"]), [-3., 5.])
        self.assertEqual(set(data.obs.donor_source), {"HPAP"})
        self.assertNotIn("phenotype", data.obs)

    def test_duplicate_object_ids_rejected(self):
        self.frame.loc[1, "Object ID"] = "a"
        with self.assertRaisesRegex(ValueError, "duplicate"):
            analysis.convert_frame(self.frame, ["INS"], self.row, self.meta)

    def test_wrong_source_image_rejected(self):
        self.frame.loc[1, "Image"] = "wrong-donor.tif"
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            analysis.convert_frame(self.frame, ["INS"], self.row, self.meta)

    def test_missing_source_image_rejected(self):
        self.frame.loc[1, "Image"] = None
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            analysis.convert_frame(self.frame, ["INS"], self.row, self.meta)

    def test_missing_or_invalid_marker_is_never_imputed(self):
        with self.assertRaisesRegex(ValueError, "Missing measurements"):
            analysis.convert_frame(self.frame, ["CD45"], self.row, self.meta)
        self.frame.loc[1, "Cell: INS: Mean"] = np.nan
        with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
            analysis.convert_frame(self.frame, ["INS"], self.row, self.meta)

    def test_identity_hold_blocks_before_file_access(self):
        with self.assertRaisesRegex(ValueError, "identity_hold"):
            analysis.import_case(Path("/not-used"), dict(self.row, inclusion_status="identity_hold"), [], {}, MODULE, MODULE)

    def test_archival_mode_retains_unmeasured_objects_without_imputation(self):
        self.frame.loc[1, ["Cell: INS: Mean", "Cell: GCG: Mean"]] = np.nan
        with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
            analysis.convert_frame(self.frame, ["INS", "GCG"], self.row, self.meta)
        data = analysis.convert_frame(self.frame, ["INS", "GCG"], self.row, self.meta, retain_unmeasured=True)
        self.assertEqual(list(data.obs.index), ["78:a", "78:b"])
        self.assertEqual(list(data.obs.marker_measurements_complete), [True, False])
        self.assertTrue(np.isnan(data.X[1]).all())
        self.assertTrue(np.isnan(data.layers["log"][1]).all())
        self.assertNotIn("phenotype", data.obs)

    def test_archival_mode_still_rejects_partial_missingness_and_invalid_values(self):
        for values in [[np.nan, 1.], [np.inf, np.nan], [-1., 1.]]:
            self.frame.loc[1, ["Cell: INS: Mean", "Cell: GCG: Mean"]] = values
            with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                analysis.convert_frame(self.frame, ["INS", "GCG"], self.row, self.meta, retain_unmeasured=True)
        self.frame[["Cell: INS: Mean", "Cell: GCG: Mean"]] = np.nan
        with self.assertRaisesRegex(ValueError, "No cells"):
            analysis.convert_frame(self.frame, ["INS", "GCG"], self.row, self.meta, retain_unmeasured=True)

    def test_numeric_spreadsheet_id_and_hpap_provenance(self):
        self.assertEqual(analysis.case_id("0078.0"), "78")
        for bad in ["nan", "78.2", "-1", "HPAP078", ""]:
            with self.assertRaises(ValueError):
                analysis.case_id(bad)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "key.csv"
            pd.DataFrame({"Case ID": [78.0], "Age": [18.], "Gender": ["F"], "Donor Status": ["T1D"]}).to_csv(path, index=False)
            result = analysis.donor_metadata(path, [{"case_id": "78", "donor_source": "HPAP"}])
            self.assertEqual(result["78"]["donor_source"], "HPAP")

    def test_complete_import_resume_and_truncation_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            exports = root / "exports/0078"
            exports.mkdir(parents=True)
            logs = root / "logs/0078"
            logs.mkdir(parents=True)
            (logs / "export.done.json").write_text("{}")
            cells = exports / "cells.tsv"
            self.frame.loc[0, "Centroid X µm"] = 3.0000000000000004
            self.frame.to_csv(cells, sep="\t", index=False)
            manifest = {"case_id": "78", "source_filename": "hp.tif", "cell_count": 2,
                        "coordinate_units": "micrometres", "pixel_size_um_x": 0.3774,
                        "pixel_size_um_y": 0.3774, "cells_sha256": analysis.sha256(cells),
                        "numeric_encoding": "java_double_to_string"}
            (exports / "export_manifest.json").write_text(json.dumps(dict(manifest, numeric_encoding="automatic_rounding")))
            with self.assertRaisesRegex(ValueError, "Full-precision"):
                analysis.import_case(root, self.row, ["INS", "GCG"], self.meta, MODULE, MODULE)
            (exports / "export_manifest.json").write_text(json.dumps(manifest))
            reports = root / "reports/image_conversion"
            reports.mkdir(parents=True)
            (reports / "input_inventory.json").write_text("[]")
            (reports / "hp.json").write_text(json.dumps({"source": "/source/hp.tif", "status": "validated"}))
            result = analysis.import_case(root, self.row, ["INS", "GCG"], self.meta, MODULE, MODULE)
            self.assertEqual(result["cells"], 2)
            imported = analysis.ad.read_h5ad(result["output"])
            self.assertEqual(imported.obs["X_centroid"].iloc[0], 3.0000000000000004)
            self.assertEqual(result, analysis.import_case(root, self.row, ["INS", "GCG"], self.meta, MODULE, MODULE))
            with self.assertRaisesRegex(ValueError, "Input or output changed"):
                analysis.import_case(root, self.row, ["GCG", "INS"], self.meta, MODULE, MODULE)
            with self.assertRaisesRegex(ValueError, "Input or output changed"):
                analysis.import_case(root, self.row, ["INS", "GCG"], self.meta, MODULE, MODULE, retain_unmeasured=True)
            self.frame.iloc[:1].to_csv(cells, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "hash differs"):
                analysis.import_case(root, self.row, ["INS", "GCG"], self.meta, MODULE, MODULE)

    def test_audit_detects_changed_measurements_and_object_order(self):
        row = dict(self.row, donor_source="HPAP")
        data = analysis.convert_frame(self.frame, ["INS", "GCG"], row, self.meta)
        conversion = {"shape_cyx": [2, 100, 100]}
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "cells.tsv"
            self.frame.to_csv(source, sep="\t", index=False)
            audit.validate_cell_roundtrip(data, source, ["INS", "GCG"], row, conversion)
            altered = data.copy()
            altered.X[0, 0] += 1
            with self.assertRaisesRegex(ValueError, "Raw H5AD values differ"):
                audit.validate_cell_roundtrip(altered, source, ["INS", "GCG"], row, conversion)
            altered = data.copy()
            altered.obs["Object ID"] = ["b", "a"]
            with self.assertRaisesRegex(ValueError, "identities/order differ"):
                audit.validate_cell_roundtrip(altered, source, ["INS", "GCG"], row, conversion)


class RegionalAuditTests(unittest.TestCase):
    def setUp(self):
        self.annotations = pd.DataFrame({"Name": ["Islet_1", "Islet_2"],
                                         "Classification": ["Islet", "Islet"],
                                         "Islet: cells total": [2, 3]})
        self.composition = pd.DataFrame({"name": ["Islet_2", "Islet_1"], "cells_total": [3, 2]})
        self.markers = pd.DataFrame({"region": ["Islet_1_core", "Islet_2_core"],
                                     "marker": ["CD45", "CD45"], "region_type": ["islet_core", "islet_core"]})
        self.targets = self.markers.rename(columns={"marker": "class", "region_type": "type"})

    def validate(self):
        return audit.validate_regions(self.annotations, self.composition, self.markers, self.targets, 2)

    def test_join_uses_islet_identity_and_rejects_changed_counts(self):
        self.assertEqual(self.validate()["islets"], 2)
        self.composition.loc[0, "cells_total"] = 4
        with self.assertRaisesRegex(ValueError, "cell counts differ"):
            self.validate()

    def test_missing_or_duplicate_regional_keys_are_rejected(self):
        original = self.markers.copy()
        self.markers.loc[1, "region"] = "Islet_3_core"
        with self.assertRaisesRegex(ValueError, "core coverage differs"):
            self.validate()
        self.markers = pd.concat([original, original.iloc[:1]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "Duplicate regional"):
            self.validate()

    def test_orphan_band_is_rejected(self):
        self.targets.loc[2] = ["Islet_3_band", "CD45", "islet_band"]
        with self.assertRaisesRegex(ValueError, "unknown islet"):
            self.validate()

    def test_spatial_counts_use_calibrated_centroids_and_preserve_overlap(self):
        from shapely.geometry import Polygon, mapping
        data = analysis.ad.AnnData(X=np.zeros((3, 1)))
        data.obs["X_centroid"] = np.array([1., 2.5, 5.]) * analysis.QUANTIFICATION_PIXEL_UM
        data.obs["Y_centroid"] = np.array([1., 1., 1.]) * analysis.QUANTIFICATION_PIXEL_UM
        geometries = [Polygon([(0, 0), (3, 0), (3, 3), (0, 3)]),
                      Polygon([(2, 0), (4, 0), (4, 3), (2, 3)])]
        features = [{"type": "Feature", "geometry": mapping(polygon),
                     "properties": {"name": f"Islet_{i+1}", "classification": {"name": "Islet"}}}
                    for i, polygon in enumerate(geometries)]
        composition = pd.DataFrame({"name": ["Islet_1", "Islet_2"], "cells_total": [2, 1]})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotations.geojson.gz"
            with gzip.open(path, "wt") as stream:
                json.dump({"type": "FeatureCollection", "features": features}, stream)
            result = audit.validate_spatial_counts(data, path, composition)
            self.assertEqual(result["unique_cells_in_cores"], 2)
            self.assertEqual(result["core_cell_membership_rows"], 3)
            self.assertEqual(result["cells_in_multiple_cores"], 1)
            composition.loc[0, "cells_total"] = 3
            with self.assertRaisesRegex(ValueError, "Spatial core count differs"):
                audit.validate_spatial_counts(data, path, composition)

    def test_exact_boundary_membership_matches_qupath_polygon_roi(self):
        from shapely.geometry import Polygon, mapping
        # QuPath PolygonROI includes the square's left/top boundary, excludes
        # its right/bottom boundary, and accepts the top-left corner.
        points = [(0., 1.), (1., 0.), (2., 1.), (1., 2.), (0., 0.), (1., 1.)]
        data = analysis.ad.AnnData(X=np.zeros((len(points), 1)))
        data.obs["X_centroid"] = np.array([p[0] for p in points]) * analysis.QUANTIFICATION_PIXEL_UM
        data.obs["Y_centroid"] = np.array([p[1] for p in points]) * analysis.QUANTIFICATION_PIXEL_UM
        feature = {"type": "Feature", "geometry": mapping(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])),
                   "properties": {"name": "Islet_1", "classification": {"name": "Islet"}}}
        composition = pd.DataFrame({"name": ["Islet_1"], "cells_total": [4]})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotations.geojson.gz"
            with gzip.open(path, "wt") as stream:
                json.dump({"type": "FeatureCollection", "features": [feature]}, stream)
            result = audit.validate_spatial_counts(data, path, composition)
        self.assertEqual(result["core_cell_membership_rows"], 4)
        checks = result["exact_boundary_checks"]
        self.assertEqual(len(checks), 5)
        self.assertEqual({row["object_id"] for row in checks if row["included"]}, {"0", "1", "4"})


if __name__ == "__main__":
    unittest.main()
