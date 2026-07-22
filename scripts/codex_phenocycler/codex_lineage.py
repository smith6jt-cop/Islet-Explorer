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
        "The RESTORE code's `>50` floor is a property of the code, not the panel, and in these "
        "native 16-bit units it does NOT sit between background and signal: 99.2% of all cells "
        "exceed it on ECAD (median 7632), 96.1% on panCK, 80.9% on CD45. It is effectively a "
        "no-op for the abundant gates and only bites the dim markers (SMA 28%, B3TUBB 44%). "
        "Treat it as a zero-drop, not a background floor."
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
