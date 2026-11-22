import sys
from pathlib import Path

import pandas as pd

UNIFORM = Path(__file__).resolve().parents[1] / "uniform"
sys.path.insert(0, str(UNIFORM))
from recluster_residuals import STRUCTURES
from resolve_tissue_calls import TYPES, UNDETECTED, resolve


def test_approved_tissue_endpoints_and_independent_evidence():
    registry = pd.read_csv(
        UNIFORM / "typing/registry/codex-uniform-32__pancreas.csv"
    ).set_index("marker")
    markers = registry.index[registry.role.isin(["defining", "shared"])]
    cases = [
        ("structure_marker_only", ["CD31"], "", "Stromal"),
        ("shared_marker_only", ["VIM"], "", "Excluded"),
        ("no_marker_above_background", [], "", UNDETECTED),
        ("", [], "", UNDETECTED),
        ("", ["CD45"], "", "Excluded"),
        ("", ["CD45", "CD3e", "CD8a"], "CD8a Tcell", "CD8a Tcell"),
        ("", ["CD20", "CD3e"], "", "Excluded"),
        ("endocrine_marker_only", ["INS"], "Beta cell", "Excluded"),
        ("", ["PGP9.5"], "", "Excluded"),
    ]
    rows = []
    for status, positive, candidate, _ in cases:
        row = {
            "start_class": "Tissue",
            "status": status,
            "n_positive": len(positive),
            "candidate": candidate,
            "named_vote_pass": bool(candidate),
            "runner_up": candidate,
            "compatible_reference_type": "",
        }
        for m in markers:
            present = m in positive
            for prefix in ("pos__", "canary__", "fold__", "covint__"):
                row[prefix + m] = float(present)
        for name, field, _ in STRUCTURES.values():
            row["in_" + name], row["crossed_frac_" + field] = False, 0.0
        rows.append(row)
    calls = resolve(pd.DataFrame(rows), registry)
    assert calls.label.tolist() == [case[-1] for case in cases]
    assert "Tissue" not in TYPES
    assert calls.loc[4, "rule"] == "no_supported_specific_lineage"
    assert calls.loc[6, "rule"] == "competing_supported_types"
    assert calls.loc[7, "rule"] == "endocrine_signal_outside_object"
