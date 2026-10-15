"""Small, deterministic checks for the bounded phenotype replay audit."""

import unittest

import numpy as np

from scripts.audit_phenotype_rule_replay import (
    agreement,
    apply_comparison_aliases,
    apply_tcell_correction,
    select_stratified_indices,
)


class PhenotypeRuleReplayTest(unittest.TestCase):
    def test_stratified_sample_covers_rare_and_common_groups_reproducibly(self):
        phenotype_codes = np.array([0, 0, 0, 0, 0, 1, 1, 1])
        donor_codes = np.array([0, 0, 0, 0, 1, 0, 0, 1])
        args = (phenotype_codes, donor_codes, ["A", "B"], ["D1", "D2"], 1, 42)
        indices, weights, strata = select_stratified_indices(*args)
        repeated, repeated_weights, repeated_strata = select_stratified_indices(*args)

        np.testing.assert_array_equal(indices, repeated)
        np.testing.assert_array_equal(weights, repeated_weights)
        self.assertEqual(strata, repeated_strata)
        self.assertEqual(len(strata), 4)
        self.assertEqual(len(indices), 4)
        self.assertEqual(sum(entry["population_n"] for entry in strata), 8)
        self.assertEqual(weights.sum(), 8)
        self.assertEqual({(int(phenotype_codes[i]), int(donor_codes[i])) for i in indices},
                         {(0, 0), (0, 1), (1, 0), (1, 1)})

    def test_alias_and_tcell_rule_do_not_reclassify_non_t_cells(self):
        markers = ["CD8a", "CD4"]
        x = np.array([[0.7, 0.2], [0.1, 0.8], [0.8, 0.9], [0.2, 0.1], [0.9, 0.1]])
        raw = np.array(["CD4 T cell", "CD8 T cell", "T cell", "T cell", "Ductal"])
        aliased = apply_comparison_aliases(raw)
        corrected = apply_tcell_correction(aliased, x, markers)
        self.assertEqual(corrected.tolist(),
                         ["CD8a Tcell", "CD4 Tcell", "T cell", "T cell", "Ductal"])
        self.assertEqual(aliased.tolist(),
                         ["CD4 Tcell", "CD8a Tcell", "T cell", "T cell", "Ductal"])
        scores = agreement(np.array(["A", "A", "B"]), np.array(["A", "B", "B"]),
                           np.array([3.0, 1.0, 6.0]))
        self.assertAlmostEqual(scores["sample_agreement"], 2 / 3)
        self.assertAlmostEqual(scores["population_weighted_agreement"], 0.9)


if __name__ == "__main__":
    unittest.main()
