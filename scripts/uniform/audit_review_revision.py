"""Characterize unsupported destinations and test CD45 recovery on held-out, named controls."""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "typing"))
import registry as RG
from review_revision import RETIRED, reconsider

# One donor per status, distinct from the index donor. Existing named calls
# provide concordance controls only; they are not independent ground truth.
RULE_CHECK = ("ND_6232", "Aab_6090", "T1D_6551")
HOLDOUT = ("ND_6479", "Aab_6558", "T1D_6550")


def audit_cohort(typing, dest):
    reg = RG.load()
    levels = pd.read_csv(typing / "levels.csv")
    residual, controls, sensitivity = [], [], []
    (dest / "remaining").mkdir(exist_ok=True)
    for path in sorted((dest / "audit").glob("*.parquet")):
        tag = path.stem
        audit = pd.read_parquet(path)
        d = pd.read_parquet(typing / "identity" / path.name).set_index("object_id")
        left = d.loc[audit.index[audit.disposition.eq("review_required")]]
        lv = levels[levels.image.eq(tag)].set_index("marker")
        markers = [m for m in lv.index if "pos__" + m in d]
        cols = [reg.column(m, lv.at[m, "measurement"]) for m in markers]
        raw = (
            pd.read_parquet(typing / "redsea" / path.name, columns=["object_id"] + cols)
            .set_index("object_id")
            .loc[left.index]
        )
        missing = raw.isna().any(axis=1)
        broad = left[[f"pos__{m}" for m in ("CD3e", "CD20", "CD68", "CD163")]].any(
            axis=1
        )
        subtype = left[["pos__CD4", "pos__CD8a"]].any(axis=1)
        groups = np.select(
            [left.label.eq("Unclassified"), broad, subtype],
            [
                left.status,
                "lineage marker present; candidate requirements failed",
                "CD4 or CD8a without a positive broad immune marker",
            ],
            "no positive broad immune marker or T subtype",
        )
        tally = pd.DataFrame(
            {"old_label": left.label, "reason": groups, "missing_measurement": missing}
        )
        residual.append(
            tally.value_counts().rename("cells").reset_index().assign(image=tag)
        )
        cov = pd.read_parquet(typing / "coverage" / path.name)
        interior = cov.set_index("object_id")
        trace = audit.loc[left.index].join(tally.drop(columns="old_label"))
        evidence_markers = (
            "CD45",
            "CD3e",
            "CD4",
            "CD8a",
            "CD20",
            "CD68",
            "CD163",
            "LGALS3",
            "INS",
            "GCG",
            "SST",
            "CHGA",
            "SMA",
            "panCK",
            "ECAD",
            "CK19",
            "VIM",
            "CD34",
            "PDPN",
            "CD31",
            "B3TUBB",
        )
        for m in evidence_markers:
            trace = trace.join(left[[p + m for p in ("pos__", "canary__", "fold__")]])
            trace["covint__" + m] = interior["covint__" + m].reindex(trace.index)
        trace.to_parquet(dest / "remaining" / path.name)
        # Test the proposed shortcut on both unresolved cells and previously
        # named non-T cells with the same absent broad immune markers.
        no_broad = ~d[["pos__" + m for m in ("CD3e", "CD20", "CD68", "CD163")]].any(
            axis=1
        )
        for marker in ("CD4", "CD8a"):
            signal = (
                no_broad
                & d["canary__" + marker].eq(1)
                & interior["covint__" + marker].reindex(d.index).gt(0.5)
            )
            for population, pool in [
                ("remaining", d.index.isin(left.index)),
                (
                    "named non-T controls",
                    ~d.label.isin(
                        [*RETIRED, "Tissue", "CD4 Tcell", "CD8a Tcell", "Excluded"]
                    ),
                ),
            ]:
                counts = d.loc[signal & pool].groupby("label").size()
                sensitivity.extend(
                    {
                        "image": tag,
                        "marker": marker,
                        "population": population,
                        "phenotype": label,
                        "cells": int(n),
                    }
                    for label, n in counts.items()
                )
        if tag not in HOLDOUT + RULE_CHECK:
            continue
        ctx = pd.read_parquet(typing / "context" / path.name)
        sample = d[
            d.label.isin(
                [
                    "CD8a Tcell",
                    "CD4 Tcell",
                    "B cell",
                    "Macrophage_LGALS3-",
                    "Macrophage_LGALS3+",
                ]
            )
        ].copy()
        truth = sample.label.copy()
        # Known CD45+ cells: simulate loss of the pan signal and test how often
        # the rescue can still recover the correct label using the selected lineage.
        sample["label"] = RETIRED[1]
        sample["pos__CD45"] = 0
        sample["canary__CD45"] = 0
        profile = reg.column("LGALS3", lv.at["LGALS3", "measurement"])
        profile_raw = pd.read_parquet(
            typing / "redsea" / path.name, columns=["object_id", profile]
        ).set_index("object_id")
        result = reconsider(
            sample.reset_index(), cov, ctx, np.isfinite(profile_raw[profile])
        )
        test = pd.DataFrame(
            {
                "original": truth,
                "destination": result.proposal,
                "rule": result.rule,
                "outcome": np.where(
                    result.proposal.eq(""),
                    "no recovery",
                    np.where(result.proposal.eq(truth), "concordant", "discordant"),
                ),
            }
        )
        controls.append(
            test.value_counts()
            .rename("cells")
            .reset_index()
            .assign(image=tag, cohort="holdout" if tag in HOLDOUT else "rule_check")
        )
    pd.DataFrame(sensitivity).to_csv(dest / "subtype_only_sensitivity.csv", index=False)
    pd.concat(residual).to_csv(dest / "remaining_evidence.csv", index=False)
    pd.concat(controls).to_csv(dest / "heldout_cd45_controls.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--typing", type=Path, required=True)
    p.add_argument("--dest", type=Path, required=True)
    audit_cohort(**vars(p.parse_args()))
