#!/usr/bin/env python3
"""Run with the same h5py/tifffile interpreter as the inventory command."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

MODULE_PATH = Path(__file__).resolve().parents[1] / "inventory_release_data.py"
SPEC = importlib.util.spec_from_file_location("inventory_release_data", MODULE_PATH)
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


class InventoryTests(unittest.TestCase):
    def test_json_and_xml_keep_names_without_values(self):
        structured = inventory.json_field_names({"subject": "PRIVATE_CANARY", "6505": {"value": 123}})
        self.assertEqual(structured["field_names"], ["subject", "value"])
        self.assertTrue(structured["dynamic_keys_omitted"])
        self.assertNotIn("PRIVATE_CANARY", json.dumps(structured))
        schema = inventory.xml_field_names('<OME><Image Name="PRIVATE_CANARY">private text</Image></OME>')
        self.assertEqual(schema["attribute_names"], ["Name"])
        self.assertNotIn("PRIVATE_CANARY", json.dumps(schema))
        self.assertNotIn("private text", json.dumps(schema))

    def test_csv_quoted_linebreak_is_one_row_and_values_stay_private(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "private-case.csv"
            path.write_text('case_id,note\n6505,"PRIVATE_CANARY\nsecond line"\n', encoding="utf-8")
            result = inventory.group_csv_schemas([path])
            self.assertEqual(result["schemas"][0]["rows"], 1)
            self.assertEqual(result["schemas"][0]["fields"], ["case_id", "note"])
            self.assertNotIn("PRIVATE_CANARY", json.dumps(result))
            self.assertNotIn("6505", json.dumps(result))
            self.assertNotIn("private-case", json.dumps(result))

    def test_h5ad_reports_schema_without_category_or_index_values(self):
        import h5py
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "source.h5ad"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                obs.attrs["_index"] = "_index"
                obs.create_dataset("_index", data=[b"PRIVATE_CANARY"])
                donor = obs.create_group("donor_status")
                donor.attrs["encoding-type"] = "categorical"
                donor.create_dataset("categories", data=[b"CATEGORY_CANARY"])
                donor.create_dataset("codes", data=[0])
            result = inventory.h5ad_schema(path)
            self.assertEqual(result["slots"]["obs"]["rows"], 1)
            self.assertEqual(result["slots"]["obs"]["field_names"], ["_index", "donor_status"])
            self.assertNotIn("PRIVATE_CANARY", json.dumps(result))
            self.assertNotIn("CATEGORY_CANARY", json.dumps(result))

    def test_tiff_metadata_and_sidecars_do_not_export_embedded_text(self):
        import numpy as np
        import tifffile
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tifffile.imwrite(root / "private-case.tiff", np.zeros((2, 2), dtype="uint8"),
                             description="PRIVATE_CANARY", metadata=None)
            (root / "private-case.json").write_text('{"note":"SIDECAR_CANARY"}', encoding="utf-8")
            result = inventory.image_inventory(root)
            self.assertEqual(result["tiff_files"], 1)
            self.assertEqual(result["files"], 2)
            self.assertIn("ImageDescription", result["embedded_text_tag_names"])
            self.assertFalse(result["pixel_arrays_read"])
            self.assertNotIn("PRIVATE_CANARY", json.dumps(result))
            self.assertNotIn("SIDECAR_CANARY", json.dumps(result))
            self.assertNotIn("private-case", json.dumps(result))


if __name__ == "__main__":
    unittest.main()
