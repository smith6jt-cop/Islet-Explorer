"""Finalize marker and cluster calls using the reviewed terminal classes (v26)."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import review_project as review
from recluster_residuals import HORMONES, STRUCTURES, constraints

UNDETECTED = "Panel undected"
CONFLICT = "Conflicting markers - excluded"
REASON = "conflicting_markers"
# These are the approved v25 named types. Tissue is solely an anatomical starting class.
TYPES = {
    **{s: ("Epithelial", "Endocrine", (m,)) for s, m in HORMONES.items()},
    "Acinar cell": ("Epithelial", "Exocrine", ("ECAD", "panCK")),
    "Ductal cell": ("Epithelial", "Exocrine", ("CK19",)),
    "T cell": ("Hematopoietic", "T", ("CD3e",)),
    "CD4 Tcell": ("Hematopoietic", "T", ("CD4",)),
    "CD8a Tcell": ("Hematopoietic", "T", ("CD8a",)),
    "B cell": ("Hematopoietic", "B", ("CD20",)),
    **{
        f"Macrophage_LGALS3{s}": ("Hematopoietic", "Myeloid", ("CD68", "CD163"))
        for s in ("+", "-")
    },
    "Fibroblast": ("Mesenchymal", "Fibroblast", ("VIM", "CD34", "PDPN")),
    "Stromal": ("Mesenchymal", "Mural", ("SMA",)),
    **{
        s: ("Neural" if s == "Neural" else "Endothelial", name, (m,))
        for s, (name, _, m) in STRUCTURES.items()
    },
}
TERMINAL = set(TYPES) | {"Excluded", UNDETECTED}
STATUSES = {
    "",
    "endocrine_marker_only",
    "no_marker_above_background",
    "shared_marker_only",
    "structure_marker_only",
}


def resolve(panel, registry):
    """Resolve only v25 Tissue rows. Keep per-cell evidence and disposition distinct."""
    if not set(panel.status).issubset(STATUSES):
        raise ValueError("Unrecognized Tissue status")
    markers = list(registry.index[registry.role.isin(["shared", "defining"])])
    pos = panel[["pos__" + m for m in markers]].set_axis(markers, axis=1).astype(bool)
    if not pos.sum(axis=1).eq(panel.n_positive).all():
        raise ValueError("Registry and per-cell positive-marker counts differ")
    local = pd.DataFrame(
        {
            m: pos[m] & panel["canary__" + m].eq(1) & panel["covint__" + m].gt(0.5)
            for m in markers
        }
    )
    destination = pd.Series("Excluded", index=panel.index)
    rule = pd.Series("no_supported_specific_lineage", index=panel.index)
    named_by = pd.Series("", index=panel.index)

    def set_call(mask, label, basis):
        destination.loc[mask], rule.loc[mask] = label, basis

    # Owner decisions precede cluster transfer. Preserve the exact source status as provenance.
    plain = panel.status.eq("")
    no_signal = panel.n_positive.eq(0)
    if not no_signal[panel.status.eq("no_marker_above_background")].all():
        raise ValueError("Panel-undetected status has a positive typing marker")
    set_call(no_signal, UNDETECTED, "no_marker_above_background")
    structure_markers = ["CD31", "B3TUBB", "GAP43"]
    only_structure = pos[structure_markers].any(axis=1) & ~pos.drop(
        columns=structure_markers + list(registry.index[registry.role.eq("shared")])
    ).any(axis=1)
    as_stroma = panel.status.eq("structure_marker_only")
    shared = panel.status.eq("shared_marker_only") | (
        plain
        & ~no_signal
        & ~pos.loc[:, registry.loc[markers].role.eq("defining")].any(axis=1)
    )
    set_call(shared, "Excluded", "shared_marker_only")
    set_call(as_stroma, "Stromal", "structure_marker_only")
    outside = panel.status.eq("endocrine_marker_only")
    rule.loc[outside] = "endocrine_signal_outside_object"
    eligible = (plain | outside) & ~(no_signal | as_stroma | shared)

    labels = list(TYPES)
    compatible = pd.DataFrame(
        constraints(panel, labels), index=panel.index, columns=labels
    )
    supported = pd.DataFrame(False, index=panel.index, columns=labels)
    endocrine_local = local[list(HORMONES.values())].any(axis=1)
    for label, (_, lineage, ms) in TYPES.items():
        runner = panel.runner_up.eq(label)
        residual_vote = panel.compatible_reference_type.eq(label)
        independent = runner | residual_vote
        localized = local[list(ms)].any(axis=1)
        strong_vote = panel.candidate.eq(label) & panel.named_vote_pass.eq(True)
        ok = (localized & independent) | (runner & residual_vote) | strong_vote
        # In a single locally supported lineage, a defining pan/lineage antigen can stand on its own.
        if (
            label in ("T cell", "B cell", "Ductal cell", "Stromal")
            or lineage == "Myeloid"
        ):
            ok |= localized
        if label in ("CD4 Tcell", "CD8a Tcell"):
            ok |= localized & (pos.CD3e | pos.CD45)
        if label == "T cell":
            ok &= pos.CD3e
        if label == "Fibroblast":
            ok &= ~pos.CD45 & ~pos[
                ["CD3e", "CD4", "CD8a", "CD20", "CD68", "CD163"]
            ].any(axis=1)
        if label in STRUCTURES:
            m = ms[0]
            rivals = [
                k
                for k in markers
                if k != m and not (label == "Neural" and k == "GAP43")
            ]
            folds = (
                panel[["fold__" + k for k in rivals]]
                .set_axis(rivals, axis=1)
                .where(pos[rivals], -np.inf)
            )
            ok &= panel["fold__" + m].gt(folds.max(axis=1))
        if label not in HORMONES:
            ok &= ~(outside & endocrine_local)
        supported[label] = eligible & compatible[label] & ok
    # A supported T subtype takes precedence over its broad T call.
    supported.loc[supported[["CD4 Tcell", "CD8a Tcell"]].any(axis=1), "T cell"] = False
    number = supported.sum(axis=1)
    rule.loc[number.gt(1)] = "competing_supported_types"
    for label, (_, _, ms) in TYPES.items():
        hit = supported[label] & number.eq(1)
        set_call(hit, label, "compatible_marker_evidence")
        named_by.loc[hit] = (
            pos.loc[hit, list(ms)].apply(lambda row: ";".join(row.index[row]), axis=1)
            if hit.any()
            else ""
        )
    set_call(plain & only_structure & number.eq(0), "Stromal", "structure_marker_only")
    return pd.DataFrame({"label": destination, "rule": rule, "named_by": named_by})


def build(prior, typing, baseline, audit, dest, source_tables=None):
    dest.mkdir(parents=True, exist_ok=True)
    for folder in ("identity", "cell_types", "qupath", "samples"):
        (dest / folder).mkdir(exist_ok=True)
    cluster = pd.read_parquet(
        prior / "assignments.parquet", filters=[("destination", "==", "Tissue")]
    )
    groups = pd.read_parquet(
        prior / "residual_cluster_members.parquet",
        columns=["image", "object_id", "compatible_reference_type"],
    )
    cluster = cluster.merge(groups, on=["image", "object_id"], validate="one_to_one")
    votes = pd.read_parquet(dest / "named_votes.parquet")
    counted = votes.control & votes.vote.ge(0.75) & votes.margin.ge(0.20)
    # Per-type performance is measured on different donors from the target fold; include Tissue controls as negatives.
    control = votes.loc[counted].assign(correct=lambda v: v.candidate.eq(v.label))
    qualified = []
    for fold, rows in votes.loc[votes.target].groupby("fold"):
        precision = (
            control.loc[control.fold.ne(fold)]
            .groupby("candidate")
            .correct.agg(["size", "mean"])
        )
        reliable = precision.index[
            precision["size"].ge(25) & precision["mean"].ge(0.90)
        ]
        qualified.append(
            rows.assign(
                named_vote_pass=rows.candidate.isin(reliable)
                & rows.vote.ge(0.75)
                & rows.margin.ge(0.20)
            )
        )
    target_votes = pd.concat(qualified)
    endo = pd.read_parquet(audit) if audit is not None else None
    if endo is not None and (
        endo.native_containing_objects.ne(0).any()
        or endo.in_islet.any()
        or endo.endocrine_object_id.notna().any()
    ):
        raise ValueError("Endocrine membership mismatch requires reconciliation")
    registry = pd.read_csv(
        Path(__file__).parent / "typing/registry/codex-uniform-32__pancreas.csv"
    ).set_index("marker")
    release = json.loads((baseline / "uniform_release.json").read_text())
    paths = sorted((prior / "identity").glob("*.parquet"))
    if {p.stem for p in paths} != set(release["identity_sha256"]):
        raise ValueError("Donor inventory mismatch")
    tallies, audits = [], []
    for path in paths:
        tag, case = path.stem, str(int(path.stem.split("_")[-1]))
        d = pd.read_parquet(path).set_index("object_id", verify_integrity=True)
        t = d.loc[d.label.eq("Tissue")].copy()
        old_status = t.status.copy()
        c = cluster.loc[cluster.image.eq(tag)].set_index(
            "object_id", verify_integrity=True
        )
        if set(t.index[t.status.eq("")]) != set(c.index):
            raise ValueError(f"{tag}: plain Tissue and cluster delta differ")
        if endo is not None:
            audited = endo.loc[endo.image.eq(tag), "object_id"]
            if audited.duplicated().any() or set(
                t.index[t.status.eq("endocrine_marker_only")]
            ) != set(audited):
                raise ValueError(f"{tag}: native endocrine audit and identity differ")
        ctx_cols = [
            k
            for name, field, _ in STRUCTURES.values()
            for k in ("in_" + name, "crossed_frac_" + field)
        ]
        ctx = (
            pd.read_parquet(
                typing / "context" / path.name,
                columns=["object_id", "endocrine_object_class", *ctx_cols],
            )
            .set_index("object_id", verify_integrity=True)
            .loc[t.index]
        )
        if not ctx.endocrine_object_class.fillna("Tissue").eq(t.start_class).all():
            raise ValueError(f"{tag}: endocrine membership differs")
        cv = (
            pd.read_parquet(typing / "coverage" / path.name)
            .set_index("object_id", verify_integrity=True)
            .loc[t.index]
        )
        if not cv.n_interior_pixels.eq(t.n_interior_pixels).all():
            raise ValueError(f"{tag}: interior measurements differ")
        t = (
            t.join(ctx[ctx_cols])
            .join(cv.filter(regex="^covint__"))
            .join(c[["runner_up", "compatible_reference_type"]].fillna(""))
        )
        v = target_votes.loc[target_votes.image.eq(tag)].set_index(
            "object_id", verify_integrity=True
        )
        if set(t.index[t.status.isin(["", "endocrine_marker_only"])]) != set(v.index):
            raise ValueError(f"{tag}: missing named-only predictions")
        t = t.join(v[["candidate", "vote", "margin", "named_vote_pass"]])
        calls = resolve(t, registry)
        for label, ids in calls.groupby("label").groups.items():
            comp, lin = TYPES[label][:2] if label in TYPES else ("", "")
            d.loc[ids, ["label", "compartment", "lineage"]] = [label, comp, lin]
        broad = calls.index[calls.rule.eq("structure_marker_only")]
        d.loc[broad, "lineage"] = (
            ""  # broad Stromal does not imply an SMA-positive mural subtype
        )
        ids = calls.index
        d.loc[ids, "excluded_reason"] = np.where(calls.label.eq("Excluded"), REASON, "")
        d.loc[
            ids,
            [
                "compartment_status",
                "lineage_status",
                "subtype_status",
                "status",
                "confirming_markers",
                "established_label",
            ],
        ] = ""
        d.loc[ids, "confirmed"] = np.nan
        d.loc[ids, "named_by"] = calls.named_by
        d.loc[broad, "named_by"] = "structure_marker_only"
        for col in ("compartment_basis", "lineage_basis"):
            d.loc[ids, col] = calls.rule
        d.loc[ids, "tags"] = (
            d.loc[ids, "tags"] + ";tissue_review_" + calls.rule
        ).str.strip(";")
        if set(d.label) - TERMINAL:
            raise ValueError(
                f"{tag}: unapproved terminal labels {set(d.label) - TERMINAL}"
            )
        d.reset_index().to_parquet(dest / "identity" / path.name, index=False)
        tab = pd.read_parquet(
            (source_tables or baseline / "cell_types") / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        if set(ids) != set(tab.index[tab.phenotype.eq("Tissue")]):
            raise ValueError(f"{tag}: assembly and Tissue inventory differ")
        for col in (
            "phenotype",
            "compartment",
            "lineage",
            "status",
            "named_by",
            "tags",
        ):
            tab.loc[ids, col] = d.loc[ids, "label" if col == "phenotype" else col]
        tab.loc[tab.phenotype.ne("Excluded")].reset_index().to_parquet(
            dest / "cell_types" / (case + ".parquet"), index=False
        )
        delta = d.loc[ids, [*review.IDENTITY_META, "label", "excluded_reason"]].copy()
        delta["confirmed"] = ""
        delta["old_label"] = "Tissue"
        delta["old_display"] = old_status.map(
            lambda s: "Tissue cell - " + review.START_STATUS[s] if s else "Tissue"
        )
        delta["display"] = calls.label.map(review.display).mask(
            calls.label.eq("Excluded"), CONFLICT
        )
        delta["tissue_resolution"] = calls.rule
        delta["named_type_candidate"] = v.candidate.reindex(ids)
        delta["named_type_vote"] = v.vote.reindex(ids)
        delta["previous_cluster_id"] = (
            c.fold.astype(str) + "/" + c.parent + "/" + c.child.astype(str)
        )
        delta["previous_cluster_score"] = c.score.reindex(ids)
        delta["previous_cluster_runner_up"] = c.runner_up.reindex(ids)
        # v25 cluster scores are evidence for the historical Tissue prediction only.
        for col in ("cluster_score", "cluster_runner_up", "cluster_id"):
            delta[col] = ""
        delta.reset_index().to_csv(
            dest / "qupath" / (tag + ".tsv"), sep="\t", index=False
        )
        (dest / "samples" / (tag + ".txt")).write_text("")
        rec = (
            delta[["old_display", "display", "tissue_resolution"]]
            .rename(
                columns={
                    "old_display": "old_label",
                    "display": "destination",
                    "tissue_resolution": "rule",
                }
            )
            .assign(image=tag)
        )
        audits.append(rec.reset_index())
        tallies.append(
            rec.groupby(["image", "old_label", "destination", "rule"])
            .size()
            .rename("cells")
            .reset_index()
        )
    displays = [review.display(k) for k in TYPES] + [
        UNDETECTED,
        CONFLICT,
        *review.EXCLUDED.values(),
        "Excluded - cell QC",
    ]
    (dest / "terminal_classes.json").write_text(
        json.dumps({"detection_classes": sorted(set(displays))}, indent=2) + "\n"
    )
    pd.concat(tallies).assign(disposition="retyped").to_csv(
        dest / "transitions.csv", index=False
    )
    pd.concat(audits, ignore_index=True).to_parquet(
        dest / "assignments.parquet", index=False
    )


def combine_project_deltas(checkpoint, dest, transition_name="tissue_transitions.csv"):
    """Compose both review steps by UUID when building directly from the v24 native project."""
    records = []
    for path in sorted((dest / "qupath").glob("*.tsv")):
        stage = pd.read_csv(
            checkpoint / "qupath" / path.name,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        ).set_index("object_id", verify_integrity=True)
        final = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).set_index(
            "object_id", verify_integrity=True
        )
        shared = final.index.intersection(stage.index)
        if not stage.loc[shared, "label"].eq(final.loc[shared, "old_label"]).all():
            raise ValueError(f"{path.stem}: delta chain has a source mismatch")
        if "old_display" not in stage:
            stage["old_display"] = stage.old_label.map(review.display)
        # A later step can change only subtype fields. Preserve earlier metadata updates
        # for shared UUIDs, and never clear unrelated native fields on other UUIDs.
        reserved = {"label", "old_label", "old_display", "display", "metadata_keys"}
        stage_fields = [c for c in stage if c not in reserved]
        final_fields = [c for c in final if c not in reserved]
        if "metadata_keys" not in stage:
            stage["metadata_keys"] = ";".join(stage_fields)
        final["metadata_keys"] = ";".join(final_fields)
        for key in stage_fields:
            if key not in final_fields:
                final.loc[shared, key] = stage.loc[shared, key]
        final.loc[shared, "metadata_keys"] = stage.loc[shared, "metadata_keys"].map(
            lambda keys: ";".join(dict.fromkeys(keys.split(";") + final_fields))
        )
        final.loc[shared, ["old_label", "old_display"]] = stage.loc[
            shared, ["old_label", "old_display"]
        ]
        merged = pd.concat([stage.drop(index=shared), final]).fillna("")
        merged.to_csv(path, sep="\t")
        frame = (
            merged.groupby(["old_display", "display"])
            .size()
            .rename("cells")
            .reset_index()
        )
        records.append(
            frame.rename(
                columns={"old_display": "old_label", "display": "destination"}
            ).assign(image=path.stem, disposition="retyped")
        )
    (dest / "transitions.csv").rename(dest / transition_name)
    pd.concat(records, ignore_index=True).to_csv(dest / "transitions.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("prior", "typing", "baseline", "audit", "dest"):
        parser.add_argument("--" + arg, type=Path, required=True)
    build(**vars(parser.parse_args()))
