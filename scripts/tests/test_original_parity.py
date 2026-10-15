"""Small synthetic preservation fixtures; never write to the project data tree."""
from contextlib import redirect_stdout
import hashlib
from io import StringIO
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))
from audit_original_parity import audit, main


def write_fixture(root, donors=(1, 1, 2), ids=("1_a", "1_b", "2_c")):
    root.mkdir()
    n = len(ids)
    obs = pd.DataFrame({"imageid": donors, "combined_islet_id": ids,
                        "measurement": np.arange(n, dtype=float),
                        "centroid_x": np.arange(n, dtype=float) * 10,
                        "count_Beta cell": np.arange(n, dtype=float),
                        "nullable": [np.nan] + [1.0] * (n - 1),
                        "dpt_pseudotime": np.arange(n, dtype=float) / n,
                        "dpt_pseudotime_quality": np.ones(n, dtype=float),
                        "leiden_0.3": [str(i % 2) for i in range(n)]}, index=ids)
    var = pd.DataFrame({"feature_type": ["protein", "protein"]}, index=["INS", "GCG"])
    matrix = np.arange(n * 2, dtype=float).reshape(n, 2)
    h5ad = ad.AnnData(matrix, obs=obs, var=var)
    h5ad.obsm["X_umap"] = matrix / 2
    h5ad.obsm["X_scVI_mean"] = matrix / 3
    h5ad.obsm["X_scVI_harmony_pt"] = matrix / 4
    h5ad.obsm["X_scVI_combined_harmony_pt"] = matrix / 5
    h5ad.obsm["phenotype_counts"] = matrix.copy()
    h5ad.obsm["phenotype_proportions"] = matrix / 10
    h5ad.obsm["morphology"] = matrix.copy()
    h5ad.obsp["connectivities"] = sparse.eye(n, format="csr")
    h5ad.varm["PCs"] = np.eye(2)
    h5ad.varp["feature_graph"] = sparse.eye(2, format="csr")
    h5ad.uns["neighbors"] = {"params": {"n_neighbors": 15, "metric": "cosine"}}
    h5ad.uns["source_annotation"] = {"gate": 10}
    h5ad.uns["iroot"] = 0
    h5ad.uns["groovy_composition_columns"] = np.array(["Case ID", "islet_id", "value"])
    h5ad.uns["groovy_composition_Case ID"] = np.asarray(donors)
    h5ad.uns["groovy_composition_islet_id"] = np.asarray(ids)
    h5ad.uns["groovy_composition_value"] = np.arange(n, dtype=float)
    for name in ("adata_ins_root.h5ad", "islet_explorer.h5ad"):
        h5ad.write_h5ad(root / name)
    pd.DataFrame({"case_id": donors, "islet_id": ids, "count": np.arange(n, dtype=float)}).to_csv(
        root / "barlow_islet_features.csv", index=False)


def hashes(root):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.iterdir() if p.is_file()}


class OriginalParityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.baseline = self.root / "baseline"
        self.candidate = self.root / "candidate"
        write_fixture(self.baseline)
        write_fixture(self.candidate)

    def tearDown(self):
        self.temp.cleanup()

    def mutate_h5ad(self, fn):
        path = self.candidate / "adata_ins_root.h5ad"
        value = ad.read_h5ad(path)
        value = fn(value)
        value.write_h5ad(path)

    def assertRejected(self, fragment):
        report = audit(self.baseline, self.candidate)
        self.assertFalse(report["passed"])
        self.assertTrue(any(fragment in c["name"] and not c["passed"] for c in report["checks"]), report)

    def test_equal_is_read_only(self):
        before = hashes(self.baseline), hashes(self.candidate)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"])
        self.assertEqual(report["derived_differences"], [])
        self.assertEqual(report["scope"], "islet_level_canonical_artifacts")
        self.assertFalse(report["source_cell_parity_checked"])
        self.assertEqual(before, (hashes(self.baseline), hashes(self.candidate)))

    def test_appended_donor_and_reordering_pass(self):
        # Add a donor while preserving every original value and graph edge.
        for path in self.candidate.iterdir():
            path.unlink()
        self.candidate.rmdir()
        write_fixture(self.candidate, (1, 1, 2, 3), ("1_a", "1_b", "2_c", "3_d"))
        for filename in ("adata_ins_root.h5ad", "islet_explorer.h5ad"):
            path = self.candidate / filename
            value = ad.read_h5ad(path)[[3, 2, 0, 1], [1, 0]].copy()
            value.uns["iroot"] = 2
            value.write_h5ad(path)
        self.assertTrue(audit(self.baseline, self.candidate)["passed"])

    def test_changed_numeric_default_exact_and_explicit_tolerance(self):
        def modify(value):
            value.X[0, 0] += 1e-8
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("/X")
        self.assertTrue(audit(self.baseline, self.candidate, atol=1e-7)["passed"])

    def test_missing_original_column(self):
        self.mutate_h5ad(lambda a: (a.obs.pop("measurement"), a)[1])
        self.assertRejected("required_columns")

    def test_changed_missingness(self):
        def modify(value):
            value.obs.loc["1_a", "nullable"] = 0
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("/nullable")

    def test_missing_original_islet(self):
        self.mutate_h5ad(lambda a: a[[0, 2]].copy())
        self.assertRejected("identity")

    def test_added_original_donor_islet(self):
        for path in self.candidate.iterdir():
            path.unlink()
        self.candidate.rmdir()
        write_fixture(self.candidate, (1, 1, 2, 1), ("1_a", "1_b", "2_c", "1_new"))
        self.assertRejected("identity")

    def test_duplicate_original_identity(self):
        def modify(value):
            value.obs.loc["1_b", "combined_islet_id"] = "1_a"
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("identity")

    def test_changed_refit_graph_edge_is_reported(self):
        def modify(value):
            value.obsp["connectivities"][0, 0] = 2
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        self.assertTrue(any("obsp/connectivities" in d["name"]
                            for d in report["derived_differences"]), report)

    def test_embedded_table_row_count_follows_appended_rows(self):
        # `<table>_n_rows` grows with appended donors, so it is checked for
        # consistency with the candidate table rather than for equality.
        for root, n in ((self.baseline, 3), (self.candidate, 3)):
            for filename in ("adata_ins_root.h5ad", "islet_explorer.h5ad"):
                path = root / filename
                value = ad.read_h5ad(path)
                value.uns["groovy_composition_n_rows"] = n
                value.write_h5ad(path)

        def append_row(value, n_rows):
            value.uns["groovy_composition_Case ID"] = np.append(value.uns["groovy_composition_Case ID"], 9)
            value.uns["groovy_composition_islet_id"] = np.append(value.uns["groovy_composition_islet_id"], "9_z")
            value.uns["groovy_composition_value"] = np.append(value.uns["groovy_composition_value"], 7.0)
            value.uns["groovy_composition_n_rows"] = n_rows
            return value

        self.mutate_h5ad(lambda a: append_row(a, 4))
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        self.assertTrue(any(c["name"].endswith("/uns/groovy_composition_n_rows/consistent")
                            and c["passed"] for c in report["checks"]))
        self.mutate_h5ad(lambda a: (a.uns.__setitem__("groovy_composition_n_rows", 3), a)[1])
        self.assertRejected("groovy_composition_n_rows/consistent")

    def test_missing_embedded_table_column(self):
        def modify(value):
            del value.uns["groovy_composition_value"]
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("comparison")

    def test_changed_fit_metadata_and_root_are_reported(self):
        def modify(value):
            value.uns["neighbors"]["params"]["n_neighbors"] = 99
            value.uns["iroot"] = 1
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        names = {d["name"] for d in report["derived_differences"]}
        self.assertIn("adata_ins_root.h5ad/uns/neighbors/params/n_neighbors", names)
        self.assertIn("adata_ins_root.h5ad/uns/iroot_islet_identity", names)

    def test_changed_fit_variable_loading_is_reported(self):
        def modify(value):
            value.varm["PCs"][0, 0] = 2
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        self.assertTrue(any("varm/PCs" in d["name"]
                            for d in report["derived_differences"]), report)

    def test_changed_refit_obs_and_umap_are_reported(self):
        def modify(value):
            value.obs.loc["1_a", "dpt_pseudotime"] = 0.9
            value.obs.loc["1_a", "leiden_0.3"] = "1"
            value.obsm["X_umap"][0, 0] = 42
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        names = {d["name"] for d in report["derived_differences"]}
        self.assertIn("adata_ins_root.h5ad/obs/dpt_pseudotime", names)
        self.assertIn("adata_ins_root.h5ad/obs/leiden_0.3", names)
        self.assertIn("adata_ins_root.h5ad/obsm/X_umap", names)

    def test_changed_donor_harmony_embeddings_are_reported(self):
        def modify(value):
            value.obsm["X_scVI_harmony_pt"][0, 0] = 42
            value.obsm["X_scVI_combined_harmony_pt"][0, 0] = 43
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        names = {d["name"] for d in report["derived_differences"]}
        self.assertIn("adata_ins_root.h5ad/obsm/X_scVI_harmony_pt", names)
        self.assertIn("adata_ins_root.h5ad/obsm/X_scVI_combined_harmony_pt", names)

    def test_missing_refit_fields_are_reported(self):
        def modify(value):
            del value.obs["dpt_pseudotime"]
            del value.obsm["X_umap"]
            del value.obsp["connectivities"]
            del value.varm["PCs"]
            del value.uns["neighbors"]
            del value.uns["iroot"]
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        missing = {d["name"] for d in report["derived_differences"]
                   if d["change"] == "missing_in_candidate"}
        self.assertEqual(missing, {
            "adata_ins_root.h5ad/obs/dpt_pseudotime",
            "adata_ins_root.h5ad/obsm/X_umap",
            "adata_ins_root.h5ad/obsp/connectivities",
            "adata_ins_root.h5ad/varm/PCs",
            "adata_ins_root.h5ad/uns/neighbors",
            "adata_ins_root.h5ad/uns/iroot",
        })

    def test_missing_nested_refit_metadata_is_reported(self):
        def modify(value):
            del value.uns["neighbors"]["params"]["metric"]
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertTrue(report["passed"], report)
        self.assertIn({"name": "adata_ins_root.h5ad/uns/neighbors/params/metric",
                       "change": "missing_in_candidate"}, report["derived_differences"])

    def test_changed_original_geometry_typing_and_scvi_fail(self):
        def modify(value):
            value.obs.loc["1_a", "centroid_x"] = 99
            value.obs.loc["1_b", "count_Beta cell"] = 99
            value.obs.loc["1_a", "dpt_pseudotime_quality"] = 0
            value.obsm["phenotype_counts"][0, 0] = 99
            value.obsm["phenotype_proportions"][0, 0] = 99
            value.obsm["X_scVI_mean"][0, 0] = 99
            value.uns["source_annotation"]["gate"] = 99
            return value
        self.mutate_h5ad(modify)
        report = audit(self.baseline, self.candidate)
        self.assertFalse(report["passed"])
        failed = {c["name"] for c in report["checks"] if not c["passed"]}
        self.assertTrue(any(name.endswith("/obs/centroid_x") for name in failed))
        self.assertTrue(any(name.endswith("/obs/count_Beta cell") for name in failed))
        self.assertTrue(any(name.endswith("/obs/dpt_pseudotime_quality") for name in failed))
        self.assertTrue(any(name.endswith("/obsm/phenotype_counts") for name in failed))
        self.assertTrue(any(name.endswith("/obsm/phenotype_proportions") for name in failed))
        self.assertTrue(any(name.endswith("/obsm/X_scVI_mean") for name in failed))
        self.assertTrue(any(name.endswith("/uns/source_annotation/gate") for name in failed))

    def test_raw_presence_is_required(self):
        path = self.baseline / "adata_ins_root.h5ad"
        baseline = ad.read_h5ad(path)
        baseline.raw = baseline.copy()
        baseline.write_h5ad(path)
        self.assertRejected("raw/present")

    def test_no_canonical_analogue_fails(self):
        (self.candidate / "adata_ins_root.h5ad").unlink()
        self.assertRejected("candidate_present")

    def test_changed_csv_and_optional_excel_original_rows(self):
        path = self.candidate / "barlow_islet_features.csv"
        frame = pd.read_csv(path)
        frame.loc[0, "count"] = 100
        frame.to_csv(path, index=False)
        self.assertRejected("barlow_islet_features.csv/count")
        for root in (self.baseline, self.candidate):
            pd.DataFrame({"Case ID": [1], "islet_id": [1], "count": [1 if root == self.baseline else 2]}).to_excel(
                root / "master_results.xlsx", index=False)
        self.assertRejected("master_results.xlsx/Sheet1/count")

    def test_report_must_be_outside_sources(self):
        before = hashes(self.baseline), hashes(self.candidate)
        with self.assertRaises(SystemExit):
            main(["--baseline-root", str(self.baseline), "--candidate-root", str(self.candidate),
                  "--report-root", str(self.candidate / "reports")])
        self.assertEqual(before, (hashes(self.baseline), hashes(self.candidate)))
        self.assertFalse((self.candidate / "reports").exists())

    def test_cli_failure_exit_and_external_report(self):
        (self.candidate / "barlow_islet_features.csv").unlink()
        output = self.root / "report"
        self.assertEqual(main(["--baseline-root", str(self.baseline), "--candidate-root", str(self.candidate),
                               "--report-root", str(output)]), 1)
        self.assertFalse(json.loads((output / "original_subset_parity.json").read_text())["passed"])

    def test_cli_identifies_derived_difference_and_unchecked_source_cells(self):
        def change_umap(value):
            value.obsm["X_umap"][0, 0] = 42
            return value
        self.mutate_h5ad(change_umap)
        output = self.root / "report"
        summary = StringIO()
        with redirect_stdout(summary):
            self.assertEqual(main(["--baseline-root", str(self.baseline),
                                   "--candidate-root", str(self.candidate),
                                   "--report-root", str(output)]), 0)
        result = json.loads(summary.getvalue())
        self.assertTrue(result["passed"])
        self.assertEqual(result["derived_difference_count"], 1)
        self.assertFalse(result["source_cell_parity_checked"])
        saved = json.loads((output / "original_subset_parity.json").read_text())
        self.assertEqual(saved["derived_difference_count"], 1)
        self.assertFalse(saved["source_cell_parity_checked"])


if __name__ == "__main__":
    unittest.main()
