"""Regression checks for the app's live composition trajectory, not research fits."""
import importlib.util
from pathlib import Path
import re
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "live_paga", ROOT / "app/shiny_app/python/paga_pseudotime.py")
trajectory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(trajectory)


class LiveTrajectoryFeatures(unittest.TestCase):
    def test_retired_policy_matches_the_app(self):
        source = (ROOT / "app/shiny_app/R/00_globals.R").read_text()
        assignment = re.search(r'PHENOTYPE_RETIRED\s*<-\s*c\((.*?)\)', source, re.S)
        self.assertIsNotNone(assignment)
        labels = set(re.findall(r'"([^"]+)"', assignment.group(1)))
        self.assertEqual(trajectory.RETIRED_PHENOTYPES, labels)

    def test_only_retained_ratios_enter_the_graph(self):
        columns = ["islet_id", "core_Beta cell", "peri_Beta cell", "core_Ductal",
                   "n_core_Beta cell", "n_endocrine_core"]
        columns += [f"{region}_{label}" for label in trajectory.RETIRED_PHENOTYPES
                    for region in ("core", "peri")]
        columns += ["core_Conflict B/T", "peri_Conflict_Endocrine/Exocrine", "core_conflict T/Myeloid"]
        self.assertEqual(trajectory.composition_features(columns),
                         ["core_Beta cell", "peri_Beta cell", "core_Ductal"])

    def test_retired_label_changes_cannot_move_islets_or_root(self):
        # Two versions differ only in retired labels, as happens when labels
        # are unavailable in one acquisition cohort. The graph and DPT must
        # be independent of those columns, including their order of magnitude.
        rng = np.random.default_rng(19)
        n = 100
        frame = pd.DataFrame({"islet_id": [f"1_Islet_{i}" for i in range(n)],
                              "case_id": 1, "donor_status": "ND",
                              "total_cells_core": 30, "n_endocrine_core": 20})
        for label in ("Alpha cell", "Beta cell", "Delta cell", "Ductal"):
            for region in ("core", "peri"):
                frame[f"{region}_{label}"] = rng.uniform(.05, .8, n)
        for label in sorted(trajectory.RETIRED_PHENOTYPES | {"Conflict B/T", "Conflict T/Myeloid"}):
            for region in ("core", "peri"):
                frame[f"{region}_{label}"] = rng.uniform(0, .1, n)
        altered = frame.copy()
        retired = [c for c in frame if c.startswith(("core_", "peri_"))
                   and c not in trajectory.composition_features(frame.columns)]
        altered.loc[:49, retired] = 0
        altered.loc[50:, retired] = rng.uniform(10, 100, (50, len(retired)))
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for i, values in enumerate((frame, altered)):
                path = Path(directory) / f"features-{i}.csv"
                values.to_csv(path, index=False)
                trajectory._FEATURES_CSV = str(path)
                trajectory._FEATURES = None
                trajectory._FIT_CACHE.clear()
                fit = trajectory.compute_pseudotime(25, 200)
                self.assertNotIn("error", fit)
                ctx = trajectory._context(25, 200)
                results.append((fit, ctx["root"], ctx["A"].obsp["connectivities"].copy()))
                self.assertIsNone(ctx["layout"], "DPT must not build the optional PAGA layout")
                pd.testing.assert_frame_equal(ctx["sub"], values, check_dtype=False)
            self.assertEqual(results[0][0]["islet_id"], results[1][0]["islet_id"])
            self.assertEqual(results[0][1], results[1][1])
            self.assertEqual((results[0][2] != results[1][2]).nnz, 0)
            np.testing.assert_array_equal(results[0][0]["pseudotime"], results[1][0]["pseudotime"])
        trajectory._FEATURES_CSV = None
        trajectory._FEATURES = None
        trajectory._FIT_CACHE.clear()


if __name__ == "__main__":
    unittest.main()
