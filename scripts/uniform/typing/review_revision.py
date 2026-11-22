"""Replay the v23 cluster and residual review using frozen donor gates and cell interiors."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

RETIRED = ("Hematopoietic CD45+", "Hematopoietic CD45-", "Unclassified")
FIBRO = ("VIM", "CD34", "PDPN")
EXOCRINE = ("ECAD", "panCK", "CK19")
HORMONES = ("INS", "GCG", "SST", "CHGA")
IMMUNE = {"T": ("CD3e", "CD4", "CD8a"), "B": ("CD20",), "Myeloid": ("CD68", "CD163")}


def displaced_lineages(rows):
    """Reapply the v23 lineage rank using native tie measurements saved before the override."""
    import registry as RG
    from identity import build_plan, placed

    if rows.empty:
        return pd.Series(index=rows.index, dtype=str)
    reg = RG.load()
    markers = [c[5:] for c in rows if c.startswith("pos__")]
    plan = build_plan(reg, markers)
    info = plan["comp"]["Hematopoietic"]
    lineages = list(info["lineages"])
    detect, count = [], []
    for lin in lineages:
        group = info["lineages"][lin]
        carried = [
            m
            for m in info["carried"]
            if not RG.split(reg.table.loc[m, "carrier_lineages"])
            or lin in RG.split(reg.table.loc[m, "carrier_lineages"])
        ]
        detect.append([plan["ix"][m] for m in group["specific"]])
        count.append(
            [
                plan["ix"][m]
                for m in dict.fromkeys(group["members"] + info["pan"] + carried)
            ]
        )
    evidence = [
        json.loads(value) if value else {} for value in rows.spatial_tie_evidence
    ]

    class RecordedPixels:
        def measure(self, requests, stage):
            measured = {}
            for i, left, right in requests:
                rec = evidence[i].get(stage)
                if rec is None:
                    continue
                groups = [[markers[j] for j in g] for g in (left, right)]
                if groups != rec["markers"]:
                    if groups[::-1] != rec["markers"]:
                        raise ValueError(
                            f"{rows.index[i]}: native lineage tie markers differ"
                        )
                    order = slice(None, None, -1)
                else:
                    order = slice(None)
                measured[i] = (
                    rec["interior_pixels"],
                    *rec["winner_pixels"][order],
                    *rec["median_scaled"][order],
                )
            return measured

    arrays = [
        rows[[prefix + m for m in markers]].to_numpy(dtype)
        for prefix, dtype in [("pos__", bool), ("canary__", bool), ("fold__", float)]
    ]
    result = placed(
        *arrays,
        detect,
        count,
        reg.const("CANARY_MARGIN"),
        reg.const("CANARY_MIN_DEPTH"),
        RecordedPixels(),
        np.arange(len(rows)),
        "lineage",
    )
    return pd.Series(
        [
            lineages[g] if g >= 0 and status != "lineage_ambiguous" else ""
            for g, status in zip(result["group"], result["status"])
        ],
        index=rows.index,
    )


def reconsider_native(identity, coverage, readings):
    """Use the selected profile reading from the classifier's aligned source table."""
    measured = np.isfinite(
        readings.set_index("object_id")[readings.attrs["columns"]["LGALS3"]]
    )
    return reconsider(identity, coverage, readings, measured)


def reconsider(identity, coverage, context, gal3_measured):
    """Return UUID-indexed decisions. Empty proposed labels explicitly require review."""
    cluster = identity.start_class.eq("EndocrineCluster")
    displaced = identity.tags.str.split(";").map(
        lambda tags: "canary_top_Hematopoietic" in tags
    )
    scope = identity.label.isin(RETIRED) | (
        cluster & displaced & identity.lineage.eq("Endocrine")
    )
    d = identity.loc[scope].set_index("object_id", verify_integrity=True)
    cv = coverage.set_index("object_id", verify_integrity=True).loc[d.index]
    ctx = context.set_index("object_id", verify_integrity=True).loc[d.index]
    if not np.array_equal(cv.n_interior_pixels, d.n_interior_pixels):
        raise ValueError("Interior count differs from v23 identity")
    a = d[
        ["label", "compartment_basis", "lineage", "status", "runner_up_compartment"]
    ].copy()
    a.columns = ["old_label", "old_basis", "old_lineage", "old_status", "runner_up"]
    for key in ("proposal", "compartment", "lineage", "named_by", "rule"):
        a[key] = ""

    def positive(m):
        return d["pos__" + m].eq(1)

    def localized(m):
        return positive(m) & d["canary__" + m].eq(1) & cv["covint__" + m].gt(0.5)

    def any_marker(markers, test):
        return pd.concat([test(m) for m in markers], axis=1).any(axis=1)

    def assign(mask, label, compartment, lineage, marker, rule):
        hit = mask & a.proposal.eq("")
        a.loc[hit, ["proposal", "compartment", "lineage", "named_by", "rule"]] = [
            label,
            compartment,
            lineage,
            marker,
            rule,
        ]

    # Object placement overwrites lineage; recover the saved pre-override competition.
    can = pd.DataFrame(
        {k: any_marker(v, lambda m: d["canary__" + m].eq(1)) for k, v in IMMUNE.items()}
    )
    plain = pd.DataFrame({k: any_marker(v, positive) for k, v in IMMUNE.items()})
    active = can.where(can.any(axis=1), plain, axis=0)
    immune = d.compartment.eq("Hematopoietic") | d.tags.str.contains(
        r"(?:^|;)canary_top_Hematopoietic(?:;|$)"
    )
    selected = d.lineage.copy()
    displaced_rows = immune & ~d.compartment.eq("Hematopoietic")
    selected.loc[displaced_rows] = displaced_lineages(d.loc[displaced_rows])
    measured = gal3_measured.loc[d.index]
    for lin, broad in [
        ("T", ("CD3e",)),
        ("B", ("CD20",)),
        ("Myeloid", ("CD68", "CD163")),
    ]:
        winner = selected.eq(lin)
        # CD45 exceptions require single-lineage corroboration and a majority
        # native interior for a lineage-defining marker above the canary level.
        admitted = immune & winner & any_marker(broad, positive)
        retired = d.label.isin(RETIRED)
        admitted &= positive("CD45") | (
            retired & active.sum(axis=1).eq(1) & any_marker(broad, localized)
        )
        if lin == "B":
            admitted &= active.sum(axis=1).eq(1)
        if lin == "T":
            cd8, cd4 = positive("CD8a"), positive("CD4")
            prefer8 = cd8 & (~cd4 | d.fold__CD8a.ge(d.fold__CD4))
            assign(
                admitted & prefer8,
                "CD8a Tcell",
                "Hematopoietic",
                lin,
                "CD8a",
                "T_subtype",
            )
            assign(
                admitted & cd4, "CD4 Tcell", "Hematopoietic", lin, "CD4", "T_subtype"
            )
            assign(
                retired & admitted & localized("CD3e"),
                "T cell",
                "Hematopoietic",
                lin,
                "CD3e",
                "CD3e_without_subtype",
            )
        elif lin == "B":
            assign(admitted, "B cell", "Hematopoietic", lin, "CD20", "B_lineage")
        else:
            for state, sign in [(True, "+"), (False, "-")]:
                assign(
                    admitted & measured & positive("LGALS3").eq(state),
                    "Macrophage_LGALS3" + sign,
                    "Hematopoietic",
                    lin,
                    "CD68;CD163",
                    "myeloid_profile",
                )
            assign(
                admitted & ~measured,
                "Macrophage",
                "Hematopoietic",
                lin,
                "CD68;CD163",
                "myeloid_without_profile",
            )
    # Limit further reconsideration to the retired outputs. A named endocrine
    # cell which fails immune naming keeps its established subtype.
    retired = d.label.isin(RETIRED)
    immune_broad = any_marker(("CD3e", "CD20", "CD68", "CD163"), positive)
    alternate = retired & ~(immune & immune_broad)
    structure = [
        ("Blood vessel", "CD31", "vessel", "Endothelium Blood Vessel", "Endothelial"),
        ("Nerve", "B3TUBB", "nerve", "Neural", "Neural"),
        ("Lymphatic", "PDPN", "lymphatic", "Endothelium Lymphatic", "Endothelial"),
    ]
    markers = [col[5:] for col in d if col.startswith("pos__")]
    for name, marker, field, label, comp in structure:
        rivals = [
            m for m in markers if m != marker and not (name == "Nerve" and m == "GAP43")
        ]
        top = pd.concat(
            [d["fold__" + m].where(positive(m), -np.inf) for m in rivals], axis=1
        ).max(axis=1)
        ok = (
            ctx["in_" + name].eq(True)
            & ctx["crossed_frac_" + field].gt(0.5)
            & positive(marker)
            & d["fold__" + marker].gt(top)
        )
        assign(alternate & ok, label, comp, name, marker, "structure_geometry")
    # Accept the original runner, or the unique independent alternative when
    # object placement hid the ranking. Require native support for the marker.
    exo = any_marker(EXOCRINE, localized)
    mural = localized("SMA")
    hidden = d.label.eq("Unclassified") & d.compartment_basis.isin(
        ["object", "established"]
    )
    epithelial = d.runner_up_compartment.eq("Epithelial") | (hidden & ~mural)
    stromal = d.runner_up_compartment.eq("Mesenchymal") | (hidden & ~exo)
    for ck19, label in [(True, "Ductal cell"), (False, "Acinar cell")]:
        for marker in ("CK19", "panCK", "ECAD"):
            assign(
                alternate
                & epithelial
                & exo
                & ~any_marker(HORMONES, positive)
                & positive("CK19").eq(ck19)
                & localized(marker),
                label,
                "Epithelial",
                "Exocrine",
                marker,
                "exocrine_alternative",
            )
    assign(
        alternate & stromal & mural,
        "Stromal",
        "Mesenchymal",
        "Mural",
        "SMA",
        "mural_alternative",
    )
    for m in FIBRO:
        assign(
            d.label.eq("Unclassified")
            & d.status.eq("shared_marker_only")
            & positive(m),
            "Fibroblast",
            "Mesenchymal",
            "Fibroblast",
            m,
            "fibroblast_pattern",
        )
    a["cd45_positive"] = positive("CD45")
    a["interior_pixels"] = d.n_interior_pixels
    a["disposition"] = np.select(
        [a.proposal.ne(""), retired], ["retyped", "review_required"], "keep"
    )
    a["x_um"], a["y_um"] = ctx.x_um, ctx.y_um
    return a


def apply_decisions(identity, audit):
    """Materialize only supported proposals, retaining evidence for owner review."""
    out = identity.copy().set_index("object_id", verify_integrity=True)
    changed = audit.loc[audit.disposition.eq("retyped")]
    ids = changed.index
    for key in ("compartment", "lineage", "named_by"):
        out.loc[ids, key] = changed[key]
    out.loc[ids, "label"] = changed.proposal
    out.loc[ids, "tags"] = (out.loc[ids, "tags"] + ";review_" + changed.rule).str.strip(
        ";"
    )
    for key in (
        "status",
        "subtype_status",
        "compartment_status",
        "lineage_status",
        "confirming_markers",
    ):
        out.loc[ids, key] = ""
    out.loc[ids, "compartment_basis"] = changed.rule
    out.loc[ids, "lineage_basis"] = changed.rule
    out.loc[ids, "confirmed"] = np.nan
    return out.reset_index()
