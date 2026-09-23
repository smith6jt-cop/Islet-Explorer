"""Small synthetic preservation fixtures; never write to the project data tree."""
import hashlib
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
                        "nullable": [np.nan] + [1.0] * (n - 1)}, index=ids)
    var = pd.DataFrame({"feature_type": ["protein", "protein"]}, index=["INS", "GCG"])
    matrix = np.arange(n * 2, dtype=float).reshape(n, 2)
    h5ad = ad.AnnData(matrix, obs=obs, var=var)
    h5ad.obsm["X_umap"] = matrix / 2
    h5ad.obsp["connectivities"] = sparse.eye(n, format="csr")
    h5ad.varm["PCs"] = np.eye(2)
    h5ad.varp["feature_graph"] = sparse.eye(2, format="csr")
    h5ad.uns["neighbors"] = {"params": {"n_neighbors": 15, "metric": "cosine"}}
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
        self.assertTrue(audit(self.baseline, self.candidate)["passed"])
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

    def test_changed_original_graph_edge(self):
        def modify(value):
            value.obsp["connectivities"][0, 0] = 2
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("obsp/connectivities")

    def test_missing_embedded_table_column(self):
        def modify(value):
            del value.uns["groovy_composition_value"]
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("comparison")

    def test_changed_fit_metadata_and_root_fail(self):
        def modify(value):
            value.uns["neighbors"]["params"]["n_neighbors"] = 99
            value.uns["iroot"] = 1
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("uns/neighbors/params/n_neighbors")
        self.assertRejected("iroot_islet_identity")

    def test_changed_variable_loading_fails(self):
        def modify(value):
            value.varm["PCs"][0, 0] = 2
            return value
        self.mutate_h5ad(modify)
        self.assertRejected("varm/PCs")

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


if __name__ == "__main__":
    unittest.main()
