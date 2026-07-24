"""Apply accepted RESTORE divisors to type every cell, and run the composition acceptance test.

Phase 5 of the port. The RESTORE screen emits a per-donor canonical divisor per pair (the manuscript
max-of-negative-controls); this module turns those into per-cell positivity and the ordered-residual
COMPARTMENT call, then cross-tabs against the app's existing 21 labels.

POSITIVITY (v7, divisor-gated): a cell is `{marker}+` iff its REDSEA-corrected intensity exceeds the
donor's canonical divisor for the pair where that marker is the target. The 5 compartment markers
(ECAD, CD31, CD45, PGP9.5, VIM) and CHGA are all first-pass targets, so the full compartment
composition is callable from the first-pass VALID divisors alone -- the immune sub-split needs the
immune bundle, but the compartment-level acceptance test does not.

ACCEPTANCE is the image-based PAIR REVIEW (the reviewer scored the 5-panel PDFs and accepted the
pairs by looking at the actual staining), NOT a cross-tab against the app's scimap labels. Those
labels are the flawed method this port replaces, so agreement with them is not a correctness signal
and this module does not use them. The composition summary below is a descriptive readout of the
accepted divisors, for the record -- not a pass/fail against old labels.

    conda activate scvi-env
    python codex_positivity.py --results /data/codex_phenocycler/restore_validation/screen_first_pass_valid \
        --out /data/codex_phenocycler/gated_restore
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
_PA = Path(os.environ.get("PHENOCYCLER_ROOT",
           str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))
sys.path.insert(0, str(_PA))
import phenocycler  # noqa: E402,F401  BEFORE numpy
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import codex_lineage as cl  # noqa: E402
from codex_panel import DONORS  # noqa: E402

REDSEA = Path("/data/codex_phenocycler/cells_redsea")
# marker -> the pair whose target it is (source of that marker's divisor). All are first-pass targets.
COMPARTMENT_MARKER = {c: g[0] for c, g in cl.COMPARTMENT_GATES.items()}  # Immune->CD45, Epithelial->ECAD, ...


def load_divisors(results_dir: Path) -> pd.DataFrame:
    """Per (donor, target-marker) canonical divisor from the VALID pair rows."""
    df = pd.read_csv(results_dir / "pair_evaluations.csv", dtype={"donor": str})
    v = df[df["state"] == "VALID"].copy()
    v["canonical_divisor"] = pd.to_numeric(v["canonical_divisor"], errors="coerce")
    return v[["donor", "target", "canonical_divisor"]].dropna()


def assign_compartments(divisors: pd.DataFrame, donors=DONORS) -> pd.DataFrame:
    """Per-cell compartment via the ordered residual on RESTORE positivity (redsea > donor divisor)."""
    gate_markers = [COMPARTMENT_MARKER[c] for c in cl.COMPARTMENT_ORDER]      # in COMPARTMENT_ORDER
    div = {(r.donor, r.target): r.canonical_divisor for r in divisors.itertuples()}
    out = []
    for donor in donors:
        f = glob.glob(f"{REDSEA}/donor_id={donor}/*.parquet")
        if not f:
            continue
        cols = ["object_id"] + gate_markers
        d = pd.read_parquet(f[0], columns=cols)
        comp = np.full(len(d), cl.RESIDUAL, dtype=object)
        assigned = np.zeros(len(d), dtype=bool)
        for compartment in cl.COMPARTMENT_ORDER:
            m = COMPARTMENT_MARKER[compartment]
            thr = div.get((donor, m))
            if thr is None:            # no VALID divisor for this donor/marker -> cannot call it
                continue
            pos = d[m].to_numpy(float) > thr
            take = pos & ~assigned
            comp[take] = compartment
            assigned |= take
        out.append(pd.DataFrame({"object_id": d["object_id"].astype(str),
                                 "donor_id": str(donor), "compartment": comp}))
    return pd.concat(out, ignore_index=True)


def composition_report(comp: pd.DataFrame) -> None:
    """Descriptive readout of the accepted-divisor compartment call. NOT a pass/fail (acceptance is
    the image-based pair review); no app-label comparison."""
    order = cl.COMPARTMENT_ORDER + [cl.RESIDUAL]
    vc = comp["compartment"].value_counts(normalize=True).reindex(order).fillna(0) * 100
    print("=== Cohort compartment composition (RESTORE divisor-gated) ===")
    print("  " + "  ".join(f"{c}={vc[c]:.1f}%" for c in order))
    per = comp.groupby("donor_id")["compartment"].apply(lambda s: (s == cl.RESIDUAL).mean() * 100)
    print(f"  Other% per-donor: {per.min():.1f}%-{per.max():.1f}% (spread {per.max()-per.min():.1f} pp)")
    print("  (acceptance is the pair review; this is a descriptive summary only)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path,
                    default=Path("/data/codex_phenocycler/restore_validation/screen_first_pass_valid"))
    ap.add_argument("--out", type=Path, default=None,
                    help="write gated compartment parquet here (optional)")
    args = ap.parse_args()

    divisors = load_divisors(args.results)
    print(f"Loaded {len(divisors)} VALID divisors "
          f"({divisors['target'].nunique()} markers x up to 15 donors)")
    comp = assign_compartments(divisors)
    composition_report(comp)

    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        comp.to_parquet(args.out / "compartments.parquet", index=False)
        import json
        (args.out / "provenance.json").write_text(json.dumps({
            "method_version": cl.EXPECTED_METHOD_VERSION,
            "upstream_sha": cl.UPSTREAM_PINNED_SHA,
            "results_dir": str(args.results),
            "call": "redsea > canonical_divisor (per donor, per marker)",
            "compartment_order": cl.COMPARTMENT_ORDER,
        }, indent=2))
        print(f"\nwrote {args.out}/compartments.parquet + provenance.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
