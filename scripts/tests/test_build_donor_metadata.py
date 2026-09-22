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
            "case_id,donor_source,inclusion_status\n6232,nPOD,selected\n78,HPAP,selected\n"
        )
        self.key.write_text(
            "Case ID,Donor Status,Age,Gender\n6232.0,ND,20,Female\n78.0,T1D,25,Male\n"
        )
        self.source.write_text(
            "case_id,age_years,BMI,HbA1c_percent,C_peptide_ng_mL,GADA_Result\n"
            "6232.0,20,22,5.2,<0.02,Positive\n78,99,99,99,99,Positive\n"
        )

    def build(self, with_source=False):
        return metadata.build_cohort_metadata(
            self.cohort, self.key, self.source if with_source else None
        ).set_index("case_id")

    def test_hpap_never_joins_numeric_npod_donor(self):
        out = self.build(with_source=True)
        hpap = out.loc["78"]
        self.assertEqual(hpap["source_donor_id"], "HPAP-078")
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
        self.assertEqual(list(out.index), ["6232", "78"])
        self.assertTrue(out["BMI"].isna().all())
        self.assertTrue((out["metadata_source"] == "donor_key_only").all())
        self.assertTrue(out["donor_key_sha256"].str.fullmatch("[a-f0-9]{64}").all())

    def test_identity_hold_is_retained_without_filtering(self):
        self.cohort.write_text(self.cohort.read_text().replace("6232,nPOD,selected", "6232,nPOD,identity_hold"))
        out = self.build()
        self.assertEqual(len(out), 2)
        self.assertEqual(out.loc["6232", "inclusion_status"], "identity_hold")

    def test_spreadsheet_ids_and_prefixed_ids_normalize(self):
        for value in [78, 78.0, "78.0", "0078", "HPAP-078"]:
            self.assertEqual(metadata.normalize_case_id(value), "78")
        for invalid in [78.5, "78.1", "78junk", None, True, "nan", ""]:
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                metadata.normalize_case_id(invalid)

    def test_duplicate_cohort_identity_is_rejected(self):
        with self.cohort.open("a") as handle:
            handle.write("0078.0,HPAP,selected\n")
        with self.assertRaisesRegex(ValueError, "duplicate normalized case_id"):
            self.build()

    def test_duplicate_donor_key_is_rejected(self):
        with self.key.open("a") as handle:
            handle.write("0078,T1D,25,Male\n")
        with self.assertRaisesRegex(ValueError, "duplicate normalized case_id"):
            self.build()

    def test_missing_key_donor_fails_without_partial_output(self):
        self.key.write_text("Case ID,Donor Status,Age,Gender\n6232,ND,20,Female\n")
        with self.assertRaisesRegex(ValueError, "missing from donor key"):
            self.build()

    def test_explicit_npod_source_must_have_exactly_one_match(self):
        self.source.write_text("case_id,age_years\n78,99\n")
        with self.assertRaisesRegex(ValueError, "exactly one nPOD row"):
            self.build(with_source=True)
        self.source.write_text("case_id,age_years\n6232,20\n6232.0,20\n")
        with self.assertRaisesRegex(ValueError, "exactly one nPOD row"):
            self.build(with_source=True)

    def test_wrong_source_assignment_and_status_are_rejected(self):
        for row in [" HPAP-078,nPOD,selected", "HDL78,HPAP,selected", "78,unknown,selected", "78,HPAP,excluded"]:
            with self.subTest(row=row):
                self.cohort.write_text("case_id,donor_source,inclusion_status\n" + row + "\n")
                with self.assertRaises(ValueError):
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
        self.assertTrue(pd.isna(result.loc["78", "npod_case_id"]))
        self.assertTrue(pd.isna(result.loc["78", "AAb_GADA"]))
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


if __name__ == "__main__":
    unittest.main()
