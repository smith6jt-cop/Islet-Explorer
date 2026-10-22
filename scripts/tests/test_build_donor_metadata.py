#!/usr/bin/env python3
"""Donor identity, clinical missingness, and original-command compatibility."""
import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import pandas as pd

MODULE_PATH = Path(__file__).resolve().parents[1] / "build_donor_metadata.py"
SPEC = importlib.util.spec_from_file_location("build_donor_metadata", MODULE_PATH)
metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metadata)


class DonorMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cohort = self.root / "cohort.csv"
        self.key = self.root / "donor_key.csv"
        self.source = self.root / "npod.csv"
        self.cohort.write_text(
            "case_id,donor_source,inclusion_status\n6232,nPOD,selected\n101,HPAP,selected\n"
        )
        self.key.write_text(
            "Case ID,Donor Status,Age,Gender\n6232.0,ND,20,Female\n101.0,T1D,25,Male\n"
        )
        self.source.write_text(
            "case_id,age_years,BMI,HbA1c_percent,C_peptide_ng_mL,GADA_Result\n"
            "6232.0,20,22,5.2,<0.02,Positive\n101,99,99,99,99,Positive\n"
        )

    def build(self, with_source=False):
        return metadata.build_cohort_metadata(
            self.cohort, self.key, self.source if with_source else None
        ).set_index("case_id")

    def test_hpap_never_joins_numeric_npod_donor(self):
        out = self.build(with_source=True)
        hpap = out.loc["101"]
        self.assertEqual(hpap["source_donor_id"], "HPAP-101")
        self.assertEqual(hpap["donor_source"], "HPAP")
        self.assertTrue(pd.isna(hpap["npod_case_id"]))
        self.assertEqual(hpap["age_years_key"], 25)
        for field in ["BMI", "age_years_npod", "n_aab", "AAb_GADA", "npod_source_sha256"]:
            self.assertTrue(pd.isna(hpap[field]), field)
        self.assertEqual(hpap["metadata_source"], "donor_key_only")
        self.assertEqual(out.loc["6232", "BMI"], 22)
        self.assertEqual(out.loc["6232", "metadata_source"], "donor_key+npod_export")
        self.assertEqual(len(out.loc["6232", "npod_source_sha256"]), 64)

    def test_unknown_antibodies_are_not_negative(self):
        out = self.build(with_source=True).loc["6232"]
        self.assertEqual(out["AAb_GADA"], True)
        self.assertTrue(pd.isna(out["AAb_IA2A"]))
        self.assertTrue(pd.isna(out["n_aab"]))
        self.assertTrue(out["C_peptide_censored"])
        self.assertEqual(out["C_peptide_ng_mL"], 0.02)

    def test_demographics_only_has_explicit_provenance(self):
        out = self.build()
        self.assertEqual(list(out.index), ["6232", "101"])
        self.assertTrue(out["BMI"].isna().all())
        self.assertTrue((out["metadata_source"] == "donor_key_only").all())
        self.assertTrue(out["donor_key_sha256"].str.fullmatch("[a-f0-9]{64}").all())

    def test_identity_hold_is_retained_without_filtering(self):
        self.cohort.write_text(self.cohort.read_text().replace("6232,nPOD,selected", "6232,nPOD,identity_hold"))
        out = self.build()
        self.assertEqual(len(out), 2)
        self.assertEqual(out.loc["6232", "inclusion_status"], "identity_hold")

    def test_spreadsheet_ids_and_prefixed_ids_normalize(self):
        for value in [101, 101.0, "101.0", "0101", "HPAP-101"]:
            self.assertEqual(metadata.normalize_case_id(value), "101")
        for invalid in [101.5, "101.1", "101junk", None, True, "nan", ""]:
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                metadata.normalize_case_id(invalid)

    def test_duplicate_cohort_identity_is_rejected(self):
        with self.cohort.open("a") as handle:
            handle.write("0101.0,HPAP,selected\n")
        with self.assertRaisesRegex(ValueError, "duplicate normalized case_id"):
            self.build()

    def test_duplicate_donor_key_is_rejected(self):
        with self.key.open("a") as handle:
            handle.write("0101,T1D,25,Male\n")
        with self.assertRaisesRegex(ValueError, "duplicate normalized case_id"):
            self.build()

    def test_missing_key_donor_fails_without_partial_output(self):
        self.key.write_text("Case ID,Donor Status,Age,Gender\n6232,ND,20,Female\n")
        with self.assertRaisesRegex(ValueError, "missing from donor key"):
            self.build()

    def test_explicit_npod_source_must_have_exactly_one_match(self):
        self.source.write_text("case_id,age_years\n101,99\n")
        with self.assertRaisesRegex(ValueError, "exactly one nPOD row"):
            self.build(with_source=True)
        self.source.write_text("case_id,age_years\n6232,20\n6232.0,20\n")
        with self.assertRaisesRegex(ValueError, "exactly one nPOD row"):
            self.build(with_source=True)

    def test_wrong_source_assignment_and_status_are_rejected(self):
        for row in [" HPAP-101,nPOD,selected", "HDL101,HPAP,selected", "101,unknown,selected", "101,HPAP,excluded"]:
            with self.subTest(row=row):
                self.cohort.write_text("case_id,donor_source,inclusion_status\n" + row + "\n")
                with self.assertRaises(ValueError):
                    self.build()

    def test_excluded_donor_is_rejected_in_every_spelling(self):
        for row in ["78,HPAP,selected", "HPAP-078,HPAP,selected", "HDL078,nPOD,selected", "0078,HPAP,identity_hold"]:
            with self.subTest(row=row):
                self.cohort.write_text("case_id,donor_source,inclusion_status\n6232,nPOD,selected\n" + row + "\n")
                with self.assertRaisesRegex(ValueError, "excluded donor"):
                    self.build()

    def test_prefixed_donor_key_cannot_cross_source_namespaces(self):
        self.key.write_text(self.key.read_text().replace("6232.0", "HPAP-6232"))
        with self.assertRaisesRegex(ValueError, "HPAP donor-key row"):
            self.build()

    def test_cli_preserves_original_output_for_explicit_cohort(self):
        with mock.patch("sys.argv", [str(MODULE_PATH), "--cohort", str(self.cohort), "--donor-key", str(self.key)]):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                metadata.main()
        self.assertEqual(error.exception.code, 2)

    def test_cli_rejects_output_alias_of_an_input(self):
        alias = self.root / "input_alias.csv"
        alias.symlink_to(self.key)
        before = self.key.read_bytes()
        argv = [str(MODULE_PATH), "--cohort", str(self.cohort), "--donor-key", str(self.key), "--output", str(alias)]
        with mock.patch("sys.argv", argv), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                metadata.main()
        self.assertEqual(error.exception.code, 2)
        self.assertEqual(self.key.read_bytes(), before)

    def test_cli_writes_only_candidate_with_hpap_nulls(self):
        output = self.root / "candidate" / "donor_metadata.csv"
        before = {path: path.read_bytes() for path in [self.cohort, self.key, self.source]}
        argv = [str(MODULE_PATH), "--cohort", str(self.cohort), "--donor-key", str(self.key),
                "--source", str(self.source), "--output", str(output)]
        with mock.patch("sys.argv", argv), contextlib.redirect_stdout(io.StringIO()):
            metadata.main()
        result = pd.read_csv(output, dtype={"case_id": str}).set_index("case_id")
        self.assertEqual(len(result), 2)
        self.assertTrue(pd.isna(result.loc["101", "npod_case_id"]))
        self.assertTrue(pd.isna(result.loc["101", "AAb_GADA"]))
        for path, contents in before.items():
            self.assertEqual(path.read_bytes(), contents)

    def test_original_cli_defaults_retain_fifteen_rows_and_schema(self):
        source = pd.DataFrame([
            {"case_id": source_id, "age_years": metadata.H5AD_AGE[app_id],
             "donor_type": "ND", "sex": "Female", "BMI": 22.0,
             "HbA1c_percent": 5.2, "C_peptide_ng_mL": 1.0,
             **{column: "Negative" for column in metadata.AAB_COLS}}
            for app_id, source_id in metadata.CODEX_DONORS.items()
        ])
        output = self.root / "legacy.csv"
        with mock.patch("sys.argv", [str(MODULE_PATH), "--output", str(output)]), \
                mock.patch.object(metadata, "find_source", return_value="original.xlsx") as find, \
                mock.patch.object(pd, "read_excel", return_value=source), \
                contextlib.redirect_stdout(io.StringIO()):
            metadata.main()
        find.assert_called_once_with(None)
        result = pd.read_csv(output, dtype={"case_id": str})
        self.assertEqual(set(result.case_id), set(metadata.CODEX_DONORS))
        self.assertEqual(len(result), 15)
        self.assertEqual(result.loc[result.case_id == "112", "npod_case_id"].iloc[0], "HDL112")
        self.assertEqual(list(result.columns), [
            "case_id", "npod_case_id", "donor_type", "age_years_npod", "age_group",
            "t1d_duration_years", "sex", "BMI", "HbA1c_percent", "C_peptide_ng_mL",
            "C_peptide_censored", *metadata.AAB_COLS.values(), "n_aab",
        ])



class UniformCohortTests(unittest.TestCase):
    """--uniform-cohort: every donor and field from the nPOD sheet (owner, 2026-10-22)."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.cohort = root / "cohort.csv"
        self.cohort.write_text("case_id,donor_status\n112,ND\n6534,T1D\n")
        self.source = root / "npod.xlsx"
        base = {"BMI": 17.0, "C_peptide_ng_mL": "1.2", "GADA_Result": "Negative", "IA_2A_Result": "Positive",
                "mIAA_Result": "Negative", "ZnT8A_Result": "Negative", "height_cm": 96, "weight_kg": 15.7,
                "BMI_percentile": None, "pancreas_weight_grams": 27.56, "admission_glucose_mg_dL": 512,
                "peak_glucose_mg_dL": 512, "RIN": 7.8, "race_ethnicity": "White", "cause_of_death": "Anoxia"}
        rows = [
            dict(base, case_id="HDL112", donor_type="No Diabetes", age_years=3.72, diabetes_hx_years=None,
                 sex="Female", HbA1c_percent=0.0, ICU_time_days="4.41 days", transit_time_minutes="960.00 minutes",
                 case_recovery_type="Organ Donor", ABO_blood_type="B+"),
            dict(base, case_id="6534", donor_type="Type 1 Diabetes", age_years=4.19, diabetes_hx_years=0.0,
                 sex="Female", HbA1c_percent=15.6, ICU_time_days="7.72 days",
                 transit_time_minutes="-1034.00 minutes", case_recovery_type="DCD Organ Donor", ABO_blood_type="O"),
        ]
        pd.DataFrame(rows).to_excel(self.source, sheet_name="data", index=False)

    def test_fields_units_and_impossible_values(self):
        out = metadata.build_uniform_metadata(self.cohort, self.source).set_index("case_id")
        self.assertEqual(out.loc["112", "npod_case_id"], "HDL112")
        self.assertEqual(out.loc["112", "age_years_npod"], 3.72)
        self.assertTrue(pd.isna(out.loc["112", "HbA1c_percent"]))        # nPOD 0.0: not measured
        self.assertTrue(pd.isna(out.loc["112", "t1d_duration_years"]))   # no diabetes: not applicable
        self.assertEqual(out.loc["6534", "t1d_duration_years"], 0.0)
        self.assertEqual(out.loc["112", "ICU_days"], 4.41)
        self.assertTrue(pd.isna(out.loc["6534", "transit_time_min"]))    # negative: impossible
        self.assertEqual(out.loc["112", "ABO_group"], "B")
        self.assertEqual(out.loc["6534", "recovery_type"], "DCD Organ Donor")
        self.assertEqual(out.loc["6534", "n_aab"], 1)

    def test_status_must_match_the_cohort(self):
        self.cohort.write_text("case_id,donor_status\n112,T1D\n6534,T1D\n")
        with self.assertRaises(SystemExit):
            metadata.build_uniform_metadata(self.cohort, self.source)


if __name__ == "__main__":
    unittest.main()
