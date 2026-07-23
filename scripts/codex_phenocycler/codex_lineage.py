"""CODEX first-pass compartment gating + RESTORE pair configuration.

Five markers plus DAPI, not the Senior 59-plex's eight, because this panel has **CD45**: the
immune union collapses to a single pan-leukocyte gate, and Senior's `CD56 & CD57` NK rescue
(which existed only to work around the absence of CD45) is unnecessary.

    CD45+   -> Immune        ECAD+  -> Epithelial     CD31+ -> Endothelial
    PGP9.5+ -> Neural        VIM+   -> Mesenchymal    else  -> Other

ORDER CONSTRAINTS (enforced by `validate_order()`, not left to comment discipline)
---------------------------------------------------------------------------------
1. **Immune before Endothelial.** CD31/PECAM1 is expressed on monocytes, neutrophils, NK and
   naive T cells, so an endothelial-first order would steal leukocytes.
2. **Neural after Epithelial.** Islet endocrine cells are PGP9.5+ and B3TUBB+, so they must be
   cleared into the epithelial bin before the neural gate runs.

Everything else is free. ECAD and CD45 may swap -- leukocytes and epithelium are cleanly
exclusive -- so clearing the 85-90% epithelial mass first to shrink the pool for every later
gate is equally valid. `COMPARTMENT_ORDER_EPITHELIAL_FIRST` is that variant.

SMA is deliberately NOT a compartment gate: SMA+ is a subset of VIM+, so it adds nothing here.
It comes back at subtyping to make the mesenchymal call positive rather than residual.

NEURAL MARKER: PGP9.5 (UCHL1) is the anchor -- the most robust pan-neuronal signal, and its
expression in islet endocrine cells is harmless because gate order removes endocrine first.
B3TUBB is equivalent and GAP43 has lower endocrine cross-reactivity but runs dimmer in mature
quiescent nerve; both are retained in `NEURAL_CONSENSUS` for a downstream consensus call.
"""

from __future__ import annotations

from codex_panel import PANEL_31, to_panel

# --- Compartments -----------------------------------------------------------------------------
COMPARTMENT_ORDER = ["Immune", "Epithelial", "Endothelial", "Neural", "Mesenchymal"]
COMPARTMENT_GATES = {
    "Immune":      ["CD45"],     # pan-leukocyte -- incl. mast and plasma cells
    "Epithelial":  ["ECAD"],     # endocrine + exocrine
    "Endothelial": ["CD31"],
    "Neural":      ["PGP9.5"],
    "Mesenchymal": ["VIM"],
}
RESIDUAL = "Other"

# The equally-valid variant: clear the epithelial mass first to shrink every later gate's pool.
COMPARTMENT_ORDER_EPITHELIAL_FIRST = ["Epithelial", "Immune", "Endothelial", "Neural", "Mesenchymal"]

SEGMENTATION_MARKER = "DAPI"          # not a phenotype gate
NEURAL_CONSENSUS = ["PGP9.5", "B3TUBB", "GAP43"]   # downstream consensus, not first pass
MESENCHYMAL_SUBTYPE = ["SMA"]         # positive muscle call at subtyping

# --- Hard ordering constraints ------------------------------------------------------------------
# (earlier, later, why). Violating one silently misassigns a whole population.
ORDER_CONSTRAINTS = [
    ("Immune", "Endothelial",
     "CD31/PECAM1 is on monocytes, neutrophils, NK and naive T cells; endothelial-first steals leukocytes"),
    ("Epithelial", "Neural",
     "islet endocrine cells are PGP9.5+/B3TUBB+ and must be cleared into Epithelial first"),
]


def validate_order(order: list[str]) -> None:
    """Raise if a candidate gate order violates a hard biological constraint."""
    missing = set(COMPARTMENT_ORDER) ^ set(order)
    if missing:
        raise ValueError(f"order must be a permutation of {COMPARTMENT_ORDER}; differs by {missing}")
    pos = {c: i for i, c in enumerate(order)}
    for earlier, later, why in ORDER_CONSTRAINTS:
        if pos[earlier] > pos[later]:
            raise ValueError(f"{earlier!r} must precede {later!r}: {why}")


# --- RESTORE pairs ------------------------------------------------------------------------------
# Directional [target, counterpart]. The counterpart's positive cells are the target's in-image
# negative control, so it must be biologically mutually exclusive with the target.
#
# Two counterparts are NOT the obvious choice, for specific reasons:
#   VIM  <- ECAD, never CD31: endothelium is vimentin-positive, so CD31 is not a valid negative.
#   PGP9.5 <- panCK, never ECAD: ECAD includes endocrine, which is PGP9.5+, so it would seed the
#            background with true neural signal. An exocrine-restricted counterpart is required.
#
# SUBSTITUTION (2026-07-22): the specified PGP9.5 counterpart was `epCAM`, which is acquired on
# only 8 of the 15 CODEX donors (the 33-channel donors 112/6521/6533/6534/6548/6549/6551 lack it).
# A counterpart present on half the cohort cannot yield a comparable cohort-wide positivity call,
# and stratifying it per donor-subset would mean two different background definitions for the same
# marker. `panCK` is the substitute: exocrine-restricted like epCAM, present on all 15 donors, and
# empirically the CLEANER negative of the two (see docs/codex_restore_pairs.md).
FIRST_PASS_PAIRS = [
    ["ECAD",   "CD31"],
    ["CD31",   "ECAD"],
    ["CD45",   "ECAD"],
    ["PGP9.5", "panCK"],   # was epCAM -- unavailable on 7/15 donors; panCK is exocrine-restricted
    ["VIM",    "ECAD"],
]

# Optional same-pass endocrine/exocrine split INSIDE the Epithelial compartment. CHGA is one
# pan-endocrine marker labelling all islet cells -- cleaner than the three hormones, and the one
# place this 31-plex is strictly better equipped than the 59-plex, which had no chromogranin.
ENDOCRINE_SPLIT_PAIR = ["CHGA", "CD31"]
ENDOCRINE_SPLIT_PARENT = "Epithelial"

# --- Input-QC arm floors (replaces the abandoned `>50` prefilter) -------------------------------
# The current-of-record method (`restore_validation.py`, METHOD_VERSION restore-pair-validation-v6)
# has dropped RESTORE's hardcoded `>50` co-positive prefilter in favour of DONOR-LOCAL,
# MARKER-SPECIFIC raw-linear floors. These are NOT divisors and NOT positivity thresholds: they
# define the provisional mutually-exclusive fitting arms, the Step 1 anchor, and the
# input-QC-retained set eligible for target-component support. Cells below both floors remain
# double-low and cannot be rescued by NNMF assignment.
#
# `estimate_marker_input_floor()` RAISES for an unregistered marker, so CODEX must declare its
# own: Senior's map has no CD45 and no PGP9.5 (its panel has neither).
#
# Assignments below are MEASURED, not inherited (all 15 donors, raw QuPath means -- see
# docs/codex_restore_pairs.md). Triangle falls below Otsu for all five gate markers, so
# max(Otsu, Triangle) reduces to plain Otsu and the second statistic buys nothing. CHGA is the
# sole exception -- Triangle EXCEEDS Otsu there, tightening 2.21% -> 1.42% high, the same
# signature that earned CD3e/CD20 their two-statistic floor in Senior.
LINEAR_OTSU_FLOOR = "linear_otsu"
MAX_LINEAR_OTSU_TRIANGLE_FLOOR = "max_linear_otsu_triangle"

MARKER_INPUT_FLOOR_METHODS = {
    "CD45":   LINEAR_OTSU_FLOOR,             # Triangle 746 << Otsu 2872
    "ECAD":   LINEAR_OTSU_FLOOR,             # Triangle 3124 << Otsu 12623
    "CD31":   LINEAR_OTSU_FLOOR,             # Triangle 603 << Otsu 3085
    "PGP9.5": LINEAR_OTSU_FLOOR,             # Triangle 543 << Otsu 5547
    "VIM":    LINEAR_OTSU_FLOOR,             # Triangle 1151 << Otsu 4460
    "panCK":  LINEAR_OTSU_FLOOR,             # Triangle 2238 < Otsu 4093 (max() moves 15.89->15.80%)
    "CHGA":   MAX_LINEAR_OTSU_TRIANGLE_FLOOR,  # Triangle 798 > Otsu 716 -- both must call it high
}

# Median donor-local floor per marker, and the min-max spread across the 15 donors. Recorded so a
# rerun that lands far outside these is caught as a regression rather than accepted silently.
# The spread is the whole reason the floors are donor-local: PGP9.5 varies 10.5x, CHGA 6.0x.
OBSERVED_FLOORS = {          # marker: (median, min, max, pct_high_at_median)
    "CD45":   (2871.7, 1812, 6653, 2.16),
    "ECAD":  (12622.9, 8384, 17170, 32.05),
    "CD31":   (3084.9, 1473, 6156, 4.42),
    "PGP9.5": (5547.0,  738, 7751, 1.70),
    "VIM":    (4459.5, 2593, 6394, 4.52),
    "panCK":  (4223.8, 3129, 7970, 15.80),
    "CHGA":    (950.9,  250, 1498, 1.42),
}


def floor_method(marker: str) -> str:
    """The arm-floor statistic for a marker. Fails closed, matching upstream."""
    m = to_panel(marker)
    if m not in MARKER_INPUT_FLOOR_METHODS:
        raise ValueError(
            f"no input-floor method is configured for marker {m!r}. Declare one here after "
            "measuring Otsu vs Triangle on this cohort -- do not inherit Senior's map blindly."
        )
    return MARKER_INPUT_FLOOR_METHODS[m]


# --- Known biases, recorded so they are not rediscovered as anomalies ---------------------------
NOTES = {
    "acinar_autofluorescence": (
        "Zymogen granules and lipofuscin make exocrine the most autofluorescent compartment, so "
        "any ECAD-high negative control OVERESTIMATES background and every ECAD-paired threshold "
        "runs conservative -- dim positives are lost, false positives are not gained. The "
        "islet-ROI mitigation does not apply: this is a whole-section compartment pass, so the "
        "negative is acinar-dominated by construction. Accept conservative first-pass cuts and "
        "refine at subtyping. For CD45 and VIM, ECAD's abundance is worth the bias."
    ),
    "prevalence": (
        "The near-zero-prevalence trap that broke KRT14/TP63/SOX2 does not bite this set -- all "
        "five gates have real prevalence in every pancreas section. Neural is sparsest "
        "(interlobular and peri-islet fibers, occasional ganglia) but reliably present, so its "
        "2-component fit will find a positive mode."
    ),
    "other_composition": (
        "Residual Other is small because CD45 sweeps up every leukocyte lineage a CD-specific "
        "union misses (mast and plasma cells now land in Immune). What remains: adipocytes (no "
        "perilipin; usually morphological), Schwann cells (no S100/SOX10 -- some caught by VIM), "
        "and vimentin-dim quiescent stellate cells."
    ),
    "intensity_floor": (
        "ABANDONED upstream. RESTORE's hardcoded `>50` co-positive prefilter survives only in the "
        "parked restore.py and as an explicit neg_q<=0 escape hatch; the current-of-record method "
        "(restore-pair-validation-v6) replaced it with donor-local marker-specific raw-linear "
        "floors -- see MARKER_INPUT_FLOOR_METHODS above. This was the right call for CODEX "
        "independently: in these native 16-bit units `>50` did NOT sit between background and "
        "signal (99.2% of cells exceed it on ECAD, median 7633; 96.1% on panCK; 80.9% on CD45), "
        "so it was a zero-drop rather than a background floor. The measured Otsu floors land at "
        "700-12,600 -- two to three orders of magnitude higher."
    ),
}


def validate_against_panel() -> None:
    """Every marker referenced here must exist in the 31-marker panel, on every donor."""
    used = {SEGMENTATION_MARKER}
    for ms in COMPARTMENT_GATES.values():
        used.update(ms)
    used.update(NEURAL_CONSENSUS)
    used.update(MESENCHYMAL_SUBTYPE)
    for tgt, cp in FIRST_PASS_PAIRS + [ENDOCRINE_SPLIT_PAIR]:
        used.update([tgt, cp])
    panel = set(PANEL_31) | {SEGMENTATION_MARKER}
    unknown = {m for m in used if to_panel(m) not in panel}
    if unknown:
        raise ValueError(
            f"markers not in the cohort-wide panel: {sorted(unknown)}. Note epCAM and KRT8-18 "
            "exist on only 8/15 donors and are deliberately excluded from cohort-wide use."
        )


def every_target_is_gated() -> None:
    """RESTORE keys a threshold on the TARGET, so every gate marker must appear at position 0."""
    targets = {t for t, _ in FIRST_PASS_PAIRS}
    gates = {m for ms in COMPARTMENT_GATES.values() for m in ms}
    ungated = gates - targets
    if ungated:
        raise ValueError(f"gate markers with no RESTORE pair (no positivity call): {sorted(ungated)}")


validate_order(COMPARTMENT_ORDER)
validate_order(COMPARTMENT_ORDER_EPITHELIAL_FIRST)
validate_against_panel()
every_target_is_gated()
