import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "uniform"))
from resolve_tcell_subtypes import BROAD, CD4, CD8, assign


def test_broad_t_cell_uses_marker_first_and_infers_only_below_both_levels():
    rows = pd.DataFrame(
        {
            "image": "i",
            "object_id": list("abcdefg"),
            "label": [BROAD] * 6 + [CD8],
            "pos__CD4": [1, 0, 1, 1, 1, 0, 0],
            "pos__CD8a": [0, 1, 1, 1, 1, 0, 1],
            "fold__CD4": [1, 0, 2, 1, 2, -1, 0],
            "fold__CD8a": [0, 1, 1, 2, 2, -1, 2],
            "cd8_vote": [0.9, 0.1, 0.9, 0.1, 0.1, 0.4, 0.1],
            "inferred_label": [CD8, CD4, CD8, CD4, CD4, CD4, CD4],
        }
    )
    calls = assign(rows)
    assert calls.destination.tolist() == [CD4, CD8, CD4, CD8, CD8, CD4]
    assert calls.rule.tolist() == ["single_positive"] * 2 + [
        "higher_of_two_positive"
    ] * 3 + ["inferred_below_both_levels"]
    assert calls.subtype_vote.isna().tolist() == [True] * 5 + [False]
    assert calls.subtype_vote.iloc[-1] == 0.6


def test_combined_review_retains_per_cell_metadata_and_original_display(tmp_path):
    from resolve_tissue_calls import combine_project_deltas

    prior = tmp_path / "prior"
    final = tmp_path / "final"
    for root in (prior, final):
        (root / "qupath").mkdir(parents=True)
    stage = [
        [
            "a",
            "Tissue",
            "Original tissue display",
            BROAD,
            "T Cell",
            "0.8",
            "",
            "cluster_score",
        ],
        [
            "b",
            "Tissue",
            "Original tissue display",
            "Stromal",
            "Stromal",
            "",
            "structure",
            "route",
        ],
    ]
    base = ["object_id", "old_label", "old_display", "label", "display"]
    pd.DataFrame(
        stage, columns=base + ["cluster_score", "route", "metadata_keys"]
    ).to_csv(prior / "qupath/image.tsv", sep="\t", index=False)
    subtype = [
        ["a", BROAD, "T Cell", CD4, "CD4 display", "inferred"],
        ["c", BROAD, "T Cell", CD8, "CD8 display", "direct"],
    ]
    path = final / "qupath/image.tsv"
    pd.DataFrame(subtype, columns=base + ["subtype_status"]).to_csv(
        path, sep="\t", index=False
    )
    (final / "transitions.csv").write_text("cells\n2\n")
    combine_project_deltas(prior, final, "t_subtype_transitions.csv")
    merged = pd.read_csv(path, sep="\t", keep_default_na=False).set_index("object_id")
    assert merged.loc[
        "a", ["old_label", "old_display", "label", "cluster_score", "subtype_status"]
    ].tolist() == ["Tissue", "Original tissue display", CD4, "0.8", "inferred"]
    assert merged.metadata_keys.to_dict() == {
        "a": "cluster_score;subtype_status",
        "b": "route",
        "c": "subtype_status",
    }
    assert merged.at["b", "route"] == "structure"
    assert merged.at["c", "old_label"] == BROAD
