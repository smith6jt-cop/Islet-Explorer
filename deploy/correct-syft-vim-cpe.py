#!/usr/bin/env python3
"""Remove one confirmed Syft 1.51.1 CRAN VIM/editor CPE collision."""

import argparse
import copy
import hashlib
import json
from pathlib import Path


CORRECTION_ID = "syft-1.51.1-cran-vim-6.2.6-editor-cpe"
INCORRECT_CPE = "cpe:2.3:a:VIM:VIM:6.2.6:*:*:*:*:*:*:*"
REQUIRED_PROPERTIES = {
    "syft:package:foundBy": "r-package-cataloger",
    "syft:package:type": "R-package",
    "syft:package:language": "R",
    "syft:package:metadataType": "r-description",
    "syft:location:0:path": "/usr/local/lib/R/site-library/VIM/DESCRIPTION",
}
SOURCES = [
    "https://cran.r-project.org/package=VIM",
    "https://github.com/statistikat/VIM",
]


def correct_document(document):
    """Return a copy with only the exact erroneous component CPE removed."""
    if document.get("bomFormat") != "CycloneDX" or not isinstance(document.get("components"), list):
        raise ValueError("Expected a CycloneDX document with a components array")
    corrected = copy.deepcopy(document)
    changes = []
    tools = document.get("metadata", {}).get("tools", {})
    syft_matches = isinstance(tools, dict) and any(
        tool.get("name") == "syft" and tool.get("version") == "1.51.1"
        for tool in tools.get("components", [])
    )
    if not syft_matches:
        return corrected, changes
    for index, component in enumerate(corrected["components"]):
        if not all(component.get(key) == value for key, value in {
            "type": "library", "name": "VIM", "version": "6.2.6",
            "purl": "pkg:cran/VIM@6.2.6", "cpe": INCORRECT_CPE,
        }.items()):
            continue
        properties = component.get("properties", [])
        if not all(
            [prop.get("value") for prop in properties if prop.get("name") == name] == [value]
            for name, value in REQUIRED_PROPERTIES.items()
        ):
            continue
        del component["cpe"]
        changes.append({"json_pointer": f"/components/{index}/cpe",
                        "removed_value": INCORRECT_CPE,
                        "component_purl": component["purl"]})
    return corrected, changes


def write_correction(source, target, record):
    paths = [path.resolve() for path in (source, target, record)]
    if len(set(paths)) != 3 or target.exists() or record.exists():
        raise ValueError("Input and new output paths must be distinct; outputs must not exist")
    original = source.read_bytes()
    document, changes = correct_document(json.loads(original))
    # A no-op preserves even the original JSON formatting.
    output = (json.dumps(document, indent=2, ensure_ascii=False) + "\n").encode() if changes else original
    evidence = {
        "schema_version": 1, "correction_id": CORRECTION_ID,
        "original_sbom_sha256": hashlib.sha256(original).hexdigest(),
        "scan_sbom_sha256": hashlib.sha256(output).hexdigest(),
        "correction_count": len(changes), "changes": changes,
        "reason": "CRAN VIM is the R missing-value visualization/imputation package; the editor CPE is an incorrect package identity.",
        "identity_sources": SOURCES,
    }
    with target.open("xb") as stream:
        stream.write(output)
    with record.open("x", encoding="utf-8") as stream:
        json.dump(evidence, stream, indent=2)
        stream.write("\n")
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("record", type=Path)
    args = parser.parse_args()
    try:
        result = write_correction(args.source, args.target, args.record)
    except (OSError, ValueError, TypeError, AttributeError) as error:
        parser.exit(1, f"SBOM identity correction failed: {error}\n")
    print(f"SBOM identity corrections: {result['correction_count']}")


if __name__ == "__main__":
    main()
