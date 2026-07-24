"""CODEX adapter for the Gate-2 reference-selection dossier (upstream reference_selection.py, v8).

Reads the two-clip-group level-2 screen, concatenates them, and writes the reference dossier
(reference_matrix.csv + frozen_reference_matrix.csv + one figure per target + manifest + SUMMARY.md),
picking ONE production reference per target. It writes NO divisor, no normalized value, and no
phenotype call -- a review artifact, not a production step.

WHY EACH PIECE EXISTS (verified against upstream reference_selection.py @ 8db37af / v8):
  * `reference_selection` binds `CANDIDATE_PAIRS`, `BASELINE_PAIRS`, `PairValidationConfig` AND
    `_TARGET_ORDER` at IMPORT (the first three from restore_validation). Patching restore_validation.*
    after import does NOT reach these already-bound copies, so `patch_reference_selection()` rebinds
    `rs.CANDIDATE_PAIRS` / `rs.BASELINE_PAIRS` / `rs._TARGET_ORDER` explicitly. Without `_TARGET_ORDER`
    the dossier writes the matrix but ZERO figures (`build_bundle` only plots targets in `_TARGET_ORDER`,
    which upstream freezes to the 7 Senior targets).
  * `build_bundle` defaults `overrides=None` -> upstream's Senior `ACCEPTED_REFERENCE_OVERRIDES`
    (`{CD11b: EpCAM, E_cadherin: CD31}`); CODEX MUST pass its own overrides.
  * The clip is irrelevant to reference_selection (it reloads no samples, uses only
    `min_reference_to_target_ratio`), so the two clip-group screens concatenate cleanly into one
    v8 `pair_evaluations.csv`. Degenerate (INPUT_QC_DEGENERATE) rows carry no valid reference
    measurement and are dropped before the matrix.

    conda activate scvi-env
    python codex_reference_selection.py     # reads screen_l2_{clipped,unclipped} -> reference_dossier_l2
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

# codex modules + phenocycler BEFORE numpy (BLAS pinning).
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
_PA = Path(os.environ.get("PHENOCYCLER_ROOT",
           str(Path(__file__).resolve().parents[2] / "external" / "Phenocycler_Analysis")))
sys.path.insert(0, str(_PA))
import phenocycler  # noqa: E402,F401  BEFORE numpy
import pandas as pd  # noqa: E402
from phenocycler import reference_selection as rs  # noqa: E402
from phenocycler import restore_validation as rv  # noqa: E402

import codex_lineage as cl  # noqa: E402
import codex_restore_validation as crv  # noqa: E402

WORKDIR = crv.WORKDIR

# One predeclared reference per target (the "baseline"). Each target appears exactly once here. The
# 3 hormone comparators (INS/GCG/SST <- panCK) are EXTRA candidates only (in crv.CODEX_PAIRS, not
# here), so for those targets baseline != candidate and the dossier can weigh the comparator.
CODEX_BASELINE_PAIRS = tuple(tuple(p) for p in (cl.FIRST_PASS_PAIRS + cl.LEVEL2_PAIRS))

# Accept-despite-undersized overrides: force the spec's reference even where a comparator passes the
# pure-proportion rule (the comparator still appears in the matrix/figures as the evidence for review).
# Each maps a target to a screened reference; recommend_references records the caveat automatically.
CODEX_REFERENCE_OVERRIDES = {
    # Epithelial-family targets -> CD31 (endothelium is the sparse-but-exclusive negative). All read
    # UNDERSIZED because epithelium >> endothelium, exactly the reviewer-accepted ECAD<-CD31 precedent
    # from level 1. The frequency rule is advisory for these structural targets; accept per the spec.
    "ECAD":   "CD31",   # level-1 precedent (reviewer-accepted, divisor 8582)
    "CK19":   "CD31",   # ductal -- epithelial family; undersized 93% like ECAD
    "panCK":  "CD31",   # epithelial confirmatory; undersized 100% like ECAD
    "LGALS3": "CD31",   # stated exception (dual macrophage+ductal -> ECAD disqualified as counterpart)
    "INS":    "GCG",    # alpha < beta in ND -> GCG undersized; the spec pairs beta as INS <- GCG
    "GCG":    "INS",    # beta-loss in T1D -> INS undersized; the spec pairs alpha as GCG <- INS
    "SST":    "INS",    # beta-loss in T1D -> INS undersized; the spec pairs delta as SST <- INS
}


def patch_reference_selection() -> None:
    """Rebind the reference_selection import-time copies to the CODEX pair web. See module docstring."""
    crv.patch_upstream()                                     # rv.CANDIDATE_PAIRS / floors / pin guard
    rs.CANDIDATE_PAIRS = tuple(crv.CODEX_PAIRS)              # REQUIRED (import-time trap)
    rs.BASELINE_PAIRS = CODEX_BASELINE_PAIRS                 # REQUIRED
    rs._TARGET_ORDER = tuple(t for t, _ in CODEX_BASELINE_PAIRS)  # REQUIRED or zero figures emitted
    # Sanity: every override target/reference is a screened pair.
    screened = set(crv.CODEX_PAIRS)
    for tgt, ref in CODEX_REFERENCE_OVERRIDES.items():
        if (tgt, ref) not in screened:
            raise ValueError(f"override {tgt} <- {ref} is not a screened CODEX pair")


def combine_screens(clipped_dir: Path, unclipped_dir: Path, out_dir: Path, *, force: bool = False) -> Path:
    """Concatenate the two clip-group pair_evaluations.csv into one v8 screen dir the dossier reads.

    Degenerate (INPUT_QC_DEGENERATE) rows carry no valid reference measurement (NaN v8 columns) and are
    dropped -- reference selection compares REFERENCE adequacy across valid fits; target positive-mode
    collapse is a positivity-step concern (codex_positivity.resolve_divisors), not a reference question.
    """
    out_dir = Path(out_dir)
    if out_dir.exists():
        if not force:
            return out_dir
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    frames = []
    for d in (Path(clipped_dir), Path(unclipped_dir)):
        pe = pd.read_csv(d / "pair_evaluations.csv")
        pe["clip_group"] = d.name
        frames.append(pe)
    combined = pd.concat(frames, ignore_index=True)
    if "reference_to_target_ratio" not in combined.columns:
        raise ValueError("screen has no reference_to_target_ratio column -- rerun at METHOD_VERSION v8")
    valid = combined[combined["reference_to_target_ratio"].notna()].copy()
    dropped = len(combined) - len(valid)
    valid.to_csv(out_dir / "pair_evaluations.csv", index=False)
    # keep one method_spec.json for provenance (reference_selection ignores the clip)
    src_spec = Path(clipped_dir) / "method_spec.json"
    if src_spec.exists():
        shutil.copy(src_spec, out_dir / "method_spec.json")
    print(f"combined {len(combined)} rows ({dropped} degenerate dropped) -> {out_dir}")
    return out_dir


def build(clipped_dir: Path, unclipped_dir: Path, output_dir: Path, *,
          work: Path = WORKDIR, force: bool = False) -> pd.DataFrame:
    patch_reference_selection()
    combined = combine_screens(clipped_dir, unclipped_dir, Path(work) / "screen_l2_combined", force=force)
    return rs.build_bundle(
        combined, Path(output_dir),
        config=rv.PairValidationConfig(),                   # min_reference_to_target_ratio = 1.0
        overrides=CODEX_REFERENCE_OVERRIDES,                # MUST pass CODEX's own, not upstream's
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clipped", type=Path, default=WORKDIR / "screen_l2_clipped")
    ap.add_argument("--unclipped", type=Path, default=WORKDIR / "screen_l2_unclipped")
    ap.add_argument("--output", type=Path, default=WORKDIR / "reference_dossier_l2")
    ap.add_argument("--work", type=Path, default=WORKDIR)
    ap.add_argument("--force-combine", action="store_true")
    args = ap.parse_args()

    matrix = build(args.clipped, args.unclipped, args.output, work=args.work, force=args.force_combine)
    cols = ["target", "reference", "ratio_median", "undersized_frac", "reference_fold_median",
            "called_frac_median", "passes_frequency_rule", "recommendation"]
    cols = [c for c in cols if c in matrix.columns]
    print("\n" + matrix[cols].to_string(index=False))
    print(f"\nwrote dossier to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
