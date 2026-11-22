"""Regression checks for the app's live composition trajectory, not research fits.

The graph inputs follow the current typing format (owner, 2026-11-18): endocrine
classes in the core only, exocrine in the peri-islet band only, macrophages by
Gal3, and no QC class.
"""
import importlib.util
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "live_paga", ROOT / "app/shiny_app/python/paga_pseudotime.py")
trajectory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(trajectory)

GLOBALS = ROOT / "app/shiny_app/R/00_globals.R"
FALLBACK = ROOT / "app/shiny_app/R/barlow_pseudotime.R"

LEGACY_MACROPHAGES = [f"{head}Macrophage_LGALS3{gal3}CD4{cd4}"
                      for head in ("CD68+CD163+", "CD68+CD163-", "CD68-CD163+")
                      for gal3 in "+-" for cd4 in "+-"]
CONFLICTS = ["Conflict B/Myeloid", "Conflict Endocrine/Exocrine", "Conflict T/Myeloid"]
# The labels of the accepted component root (twelve macrophage labels, one Exocrine class).
COMPONENT_LABELS = (["CHGA+ hormone-negative", "Beta cell", "Alpha cell", "Delta cell", "Exocrine",
                     "Hematopoietic CD45+", "Hematopoietic CD45-", "T cell", "CD8a Tcell", "CD4 Tcell",
                     "B cell", "Stromal", "Endothelium Blood Vessel", "Neural", "Endothelium Lymphatic"]
                    + LEGACY_MACROPHAGES + CONFLICTS + ["Doublet suspected", "Tissue", "Unclassified"])
# The labels the current typing code writes (two macrophage labels, Acinar / Ductal, Fibroblast).
REVISED_LABELS = (["CHGA+ hormone-negative", "Beta cell", "Alpha cell", "Delta cell", "Acinar cell",
                   "Ductal cell", "Hematopoietic CD45+", "Hematopoietic CD45-", "CD8a Tcell", "CD4 Tcell",
                   "B cell", "Stromal", "Fibroblast", "Endothelium Blood Vessel", "Neural",
                   "Endothelium Lymphatic"] + CONFLICTS
                  + ["Doublet suspected", "Macrophage_LGALS3+", "Macrophage_LGALS3-", "Tissue", "Unclassified"])
BOTH_REGIONS = ["CD8a Tcell", "CD4 Tcell", "B cell", "Stromal", "Endothelium Blood Vessel", "Neural",
                "Endothelium Lymphatic"]


def r_strings(source, name):
    """The quoted strings of an R `name <- c(...)` assignment."""
    assignment = re.search(rf'^{re.escape(name)}\s*<-\s*c\((.*?)\)', source, re.S | re.M)
    if assignment is None:
        raise AssertionError(f"{name} not found")
    return re.findall(r'"([^"]+)"', assignment.group(1))


def feature_table(labels, n=120, seed=19):
    """A feature table as assemble_uniform.py writes it: ratios and counts per label and region."""
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame({"islet_id": [f"1_Islet_{i}" for i in range(n)], "case_id": 1,
                          "donor_status": np.where(np.arange(n) % 3 == 0, "ND", "T1D"),
                          "object_class": "Islet", "total_cells_core": 40, "total_cells_peri": 60})
    endocrine = [label for label in labels if label in trajectory.ENDOCRINE_CLASSES]
    counts = {}
    for label in labels:
        for region in ("core", "peri"):
            high = 30 if label in endocrine and region == "core" else 5
            counts[f"n_{region}_{label}"] = rng.integers(0, high, n)
    counts = pd.DataFrame(counts)
    frame["n_endocrine_core"] = counts[[f"n_core_{label}" for label in endocrine]].sum(axis=1)
    denominator = frame.n_endocrine_core.clip(lower=1).to_numpy()
    ratios = {name[2:]: np.log1p(values.to_numpy() / denominator) for name, values in counts.items()}
    table = pd.concat([frame, pd.DataFrame(ratios), counts], axis=1)
    names, matrix = trajectory.phenotype_inputs(table)
    overlaps = {f"overlap__{name}__{s}": np.expm1(matrix[:, i])*denominator*rng.random(n)
                for i,name in enumerate(names) for s in trajectory.OVERLAP_STRUCTURES}
    return pd.concat([table, pd.DataFrame(overlaps)], axis=1)


def fit(frame, directory, name):
    path = Path(directory) / f"{name}.csv"
    frame.to_csv(path, index=False)
    trajectory._FEATURES_CSV = str(path)
    trajectory._FEATURES = None
    trajectory._FIT_CACHE.clear()
    result = trajectory.compute_pseudotime(25, 200)
    return result, trajectory._context(25, 200)


def reset_module():
    trajectory._FEATURES_CSV = None
    trajectory._FEATURES = None
    trajectory._FIT_CACHE.clear()


class LiveTrajectoryFeatures(unittest.TestCase):
    def tearDown(self):
        reset_module()

    def test_review_panel_results_do_not_enter_the_trajectory(self):
        for region in ("core", "peri"):
            for label in ("Panel undected", "Conflicting markers - excluded"):
                self.assertFalse(trajectory.trajectory_feature_allowed(region, label))

    def test_label_sets_match_the_app(self):
        app = GLOBALS.read_text()
        fallback = FALLBACK.read_text()
        self.assertEqual(trajectory.RETIRED_PHENOTYPES, set(r_strings(app, "PHENOTYPE_RETIRED")))
        self.assertEqual(tuple(r_strings(app, "MACROPHAGE_LABELS")), trajectory.MACROPHAGE_LABELS)
        pattern = re.search(r'^MACROPHAGE_LABEL_PATTERN\s*<-\s*"(.*)"$', app, re.M)
        self.assertEqual(pattern.group(1), trajectory.MACROPHAGE_LABEL_PATTERN)
        # names(UNCLASSIFIED_RAW_TO_KEY): the second c(...) of the setNames() call.
        family = re.search(r'UNCLASSIFIED_RAW_TO_KEY\s*<-\s*stats::setNames\(.*?\)\]\),\s*c\((.*?)\)\)',
                           app, re.S)
        self.assertEqual(trajectory.QC_CLASSES, set(re.findall(r'"([^"]+)"', family.group(1))))
        self.assertEqual(trajectory.ENDOCRINE_CLASSES,
                         set(r_strings(fallback, "TRAJECTORY_ENDOCRINE_CLASSES")))
        self.assertEqual(trajectory.EXOCRINE_CLASSES,
                         set(r_strings(fallback, "TRAJECTORY_EXOCRINE_CLASSES")))
        policy = re.search(r'^TRAJECTORY_INPUT_POLICY\s*<-\s*"(.*)"$', fallback, re.M)
        self.assertEqual(policy.group(1), trajectory.TRAJECTORY_INPUT_POLICY)

    def test_component_root_inputs_follow_the_current_typing_format(self):
        names, values = trajectory.phenotype_inputs(feature_table(COMPONENT_LABELS))
        expected = ([f"core_{label}" for label in ("CHGA+ hormone-negative", "Beta cell", "Alpha cell",
                                                   "Delta cell")] + ["peri_Exocrine"]
                    + [f"combined_{label}" for label in BOTH_REGIONS]
                    + [f"combined_{label}" for label in trajectory.MACROPHAGE_LABELS])
        self.assertEqual(names, expected)
        self.assertEqual(values.shape, (120, 14))
        self.assertTrue(np.isfinite(values).all())

    def test_revised_root_inputs_follow_the_current_typing_format(self):
        names, _ = trajectory.phenotype_inputs(feature_table(REVISED_LABELS))
        self.assertEqual(len(names), 16)
        for label in trajectory.ENDOCRINE_CLASSES:
            self.assertIn(f"core_{label}", names)
            self.assertNotIn(f"peri_{label}", names)
        for label in ("Acinar cell", "Ductal cell"):
            self.assertIn(f"peri_{label}", names)
            self.assertNotIn(f"core_{label}", names)
        for label in ("Fibroblast",) + trajectory.MACROPHAGE_LABELS:
            self.assertIn(f"combined_{label}", names)

    def test_no_qc_or_retired_class_is_an_input(self):
        labels = (sorted(trajectory.QC_CLASSES | trajectory.RETIRED_PHENOTYPES)
                  + ["Conflict B/T", "Conflict_Endocrine/Exocrine", "conflict T/Myeloid", "Conflict"])
        columns = ["islet_id", "n_endocrine_core", "n_core_Beta cell", "core_Beta cell", "peri_Beta cell",
                   "core_Exocrine", "peri_Exocrine", "core_Ductal", "peri_Ductal", "core_Macrophage"]
        columns += [f"{region}_{label}" for label in labels for region in ("core", "peri")]
        self.assertEqual(trajectory.composition_features(columns),
                         ["core_Beta cell", "peri_Exocrine", "peri_Ductal", "core_Macrophage"])
        for label in labels:
            self.assertFalse(trajectory.trajectory_feature_allowed("core", label), label)
            self.assertFalse(trajectory.trajectory_feature_allowed("peri", label), label)

    def test_macrophages_are_summed_by_gal3_before_the_log(self):
        frame = feature_table(COMPONENT_LABELS)
        names, values = trajectory.phenotype_inputs(frame)
        denominator = frame.n_endocrine_core.clip(lower=1)
        for sign in "+-":
            members = [f"n_{region}_{label}" for region in ("core","peri")
                       for label in LEGACY_MACROPHAGES if f"LGALS3{sign}" in label]
            expected = np.log1p(frame[members].sum(axis=1) / denominator)
            np.testing.assert_allclose(values[:, names.index(f"combined_Macrophage_LGALS3{sign}")], expected)
        # A class with one stored column keeps the stored ratio.
        np.testing.assert_array_equal(values[:, names.index("core_Beta cell")],
                                      frame["core_Beta cell"].to_numpy())
        revised = feature_table(REVISED_LABELS)
        names, values = trajectory.phenotype_inputs(revised)
        expected = np.log1p((revised["n_core_Macrophage_LGALS3+"] + revised["n_peri_Macrophage_LGALS3+"]) /
                            revised.n_endocrine_core.clip(lower=1))
        np.testing.assert_allclose(values[:, names.index("combined_Macrophage_LGALS3+")], expected)
        with self.assertRaisesRegex(ValueError, "needs their counts"):
            trajectory.phenotype_inputs(frame.drop(columns=[c for c in frame if c.startswith("n_peri_")]))

    def test_excluded_columns_cannot_move_objects_or_root(self):
        # Two versions differ only in what the policy leaves out: peri-islet endocrine, core exocrine,
        # QC, Conflict, and how macrophages of one Gal3 sign are divided among the older marker
        # profiles. The graph, the root and the pseudotime must be identical.
        frame = feature_table(COMPONENT_LABELS)
        altered = frame.copy()
        rng = np.random.default_rng(7)
        inputs = set(trajectory.composition_features(frame.columns))
        left_out = [c for c in frame if c.startswith(("core_", "peri_")) and c not in inputs]
        self.assertTrue({"peri_Beta cell", "core_Exocrine", "core_Unclassified", "peri_Tissue",
                         "core_Doublet suspected", "peri_Conflict T/Myeloid",
                         "core_Hematopoietic CD45+"} <= set(left_out))
        altered[left_out] = rng.uniform(10, 100, (len(frame), len(left_out)))
        altered[["n_" + c for c in left_out]] = rng.integers(50, 500, (len(frame), len(left_out)))
        for region in ("core", "peri"):
            for sign in "+-":
                members = [f"n_{region}_{label}" for label in LEGACY_MACROPHAGES
                           if f"LGALS3{sign}" in label]
                total = frame[members].sum(axis=1)
                altered[members] = 0
                altered[members[3]] = total
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for i, values in enumerate((frame, altered)):
                result, ctx = fit(values, directory, f"features-{i}")
                self.assertNotIn("error", result)
                self.assertEqual(result["n_without_inputs"], 0)
                self.assertEqual(result["input_policy"], trajectory.TRAJECTORY_INPUT_POLICY)
                self.assertEqual(len(result["features"]), 56)
                self.assertIsNone(ctx["layout"], "DPT must not build the optional PAGA layout")
                pd.testing.assert_frame_equal(ctx["sub"], values, check_dtype=False)
                results.append((result, ctx["root"], ctx["A"].obsp["connectivities"].copy()))
            self.assertEqual(results[0][0]["islet_id"], results[1][0]["islet_id"])
            self.assertEqual(results[0][0]["features"], results[1][0]["features"])
            self.assertEqual(results[0][1], results[1][1])
            self.assertEqual((results[0][2] != results[1][2]).nnz, 0)
            np.testing.assert_array_equal(results[0][0]["pseudotime"], results[1][0]["pseudotime"])

    def test_objects_without_an_input_class_get_no_pseudotime(self):
        frame = feature_table(COMPONENT_LABELS)
        names, _ = trajectory.phenotype_inputs(frame)
        sources = [c for c in trajectory.composition_features(frame.columns)]
        frame.loc[[4, 9], frame.columns.str.startswith("overlap__")] = 0
        frame.loc[[4, 9], sources] = 0.0
        frame.loc[[4, 9], ["n_" + c for c in sources]] = 0
        frame.loc[[4, 9], "n_endocrine_core"] = 0
        self.assertGreater(frame.loc[4, "n_core_Unclassified"] + frame.loc[4, "n_peri_Tissue"], 0)
        with tempfile.TemporaryDirectory() as directory:
            result, ctx = fit(frame, directory, "features")
        self.assertNotIn("error", result)
        self.assertEqual(result["n_without_inputs"], 2)
        self.assertEqual(len(result["islet_id"]), len(frame) - 2)
        self.assertFalse({"1_Islet_4", "1_Islet_9"} & set(result["islet_id"]))
        self.assertTrue(np.isfinite(result["pseudotime"]).all())

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
    def test_destiny_fallback_builds_the_same_inputs(self):
        script = r'''
            args <- commandArgs(trailingOnly = TRUE)
            env <- new.env(parent = globalenv())
            need <- c("PHENOTYPE_RETIRED", "phenotype_is_conflict", "UNCLASSIFIED_KEYS",
                      "UNCLASSIFIED_RAW_TO_KEY", "MACROPHAGE_LABELS", "MACROPHAGE_LABEL_PATTERN",
                      "MACROPHAGE_KEY_PATTERN", "is_macrophage_label", "macrophage_profile")
            for (expr in parse(args[1])) {
              if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
                  as.character(expr[[2]])[1] %in% need) eval(expr, envir = env)
            }
            sys.source(args[2], envir = env)
            d <- utils::read.csv(args[3], check.names = FALSE)
            m <- env$barlow_trajectory_inputs(d)
            utils::write.csv(as.data.frame(m, check.names = FALSE), args[4], row.names = FALSE)
        '''
        for labels in (COMPONENT_LABELS, REVISED_LABELS):
            frame = feature_table(labels)
            names, values = trajectory.trajectory_inputs(frame)
            with tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "features.csv", Path(directory) / "inputs.csv"
                frame.to_csv(source, index=False)
                done = subprocess.run(["Rscript", "-e", script, str(GLOBALS), str(FALLBACK),
                                       str(source), str(output)], capture_output=True, text=True)
                self.assertEqual(done.returncode, 0, done.stderr)
                fallback = pd.read_csv(output)
            self.assertEqual(list(fallback.columns), names)
            np.testing.assert_allclose(fallback.to_numpy(), values, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
