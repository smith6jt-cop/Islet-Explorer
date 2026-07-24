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
exclusive -- so clearing the epithelial mass first to shrink the pool for every later gate is
equally valid (order-invariant on the residual: 0.00 pp of Other). `COMPARTMENT_ORDER_EPITHELIAL_FIRST`
is that variant. (NB: a pancreas section is ~72% epithelial by the app's labels, not the "85-90%"
once stated here; and at donor-local FLOORS only ~28% reaches Epithelial -- floors are not positivity
thresholds, see NOTES["other_composition"]. RESTORE supplies the actual call.)

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

# A documented EQUIVALENCE, not a live decision: clearing the epithelial mass first is order-invariant
# on the residual, so it changes Other% by EXACTLY 0.00 pp and moves <= 0.42 pp of any single donor's
# cells (only the CD45+/ECAD+ intersection). Measured 2026-07-23 on both raw and REDSEA parquets. Do
# not spend a review round choosing between the two orders.
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

# --- Markers excluded from RESTORE in BOTH roles ------------------------------------------------
# EMPTY at level 2 (2026-07-24): CD20 is RE-INCLUDED as a B target (`CD20 <- ECAD`, in LEVEL2_PAIRS)
# per the authoritative level-2 spec. This is a deliberate CODEX-side UNWIND of the pinned upstream
# both-role exclusion. It is runtime-safe: upstream's both-role import guard runs over UPSTREAM's own
# `CANDIDATE_PAIRS` (which still excludes CD20), and `codex_restore_validation.patch_upstream()`
# rebinds `rv.CANDIDATE_PAIRS` to the CODEX set that INCLUDES `CD20 <- ECAD` -- so upstream never
# blocks the CD20 pair at runtime. The adapter still asserts upstream's `rv.RESTORE_EXCLUDED_MARKERS`
# carries CD20 (that assertion guards the SHA pin, not the CODEX policy).
#
# HISTORICAL RATIONALE for the earlier exclusion (kept for the record, now superseded): CD20 is a
# weak, donor-variable channel (donor-local prevalence 0.32%-4.85%, 15.2x swing; display-ceiling
# median ~503 counts vs ~6,800 for E-cadherin); of 42,183 CD20+ cells only 5.8% are CD45+, and 48.6%
# of the CD20-only cells are ECAD+ stray signal on epithelium (see CD20_STRAY_SIGNAL below). The
# level-2 spec accepts this: CD20 is screened as a target BUT the B call is made only INSIDE the
# Immune (CD45+) compartment at the bottom of the sub-split precedence, and CD20 is a designated
# PREVALENCE_COLLAPSE marker (borrowed/global divisor where its per-donor fit collapses).
RESTORE_EXCLUDED_MARKERS: dict[str, str] = {}

# No markers are extraction-only now that CD20 is screened. Kept for the adapter contract.
DEFERRED_EXTRACTION_MARKERS = tuple(RESTORE_EXCLUDED_MARKERS)

# The upstream commit this CODEX config is aligned to. The adapter (codex_restore_validation.py)
# asserts the consumed phenocycler package matches this, so upstream drift fails loudly rather than
# silently re-screening CD20 or reverting the both-role guard.
# 92f379a -> 55b904d (2026-07-23): audited, doc-only -- corrects the CD3e<-CD163 review caption
# (confirms CD163 is a weak-but-valid reference, reference_fold 4.4x); no pair/marker/floor/scoring
# change, so the CODEX contract is unaffected.
# 55b904d -> 58d3738 (2026-07-23): commits the v7 divisor-gated call (was the uncommitted working-tree
# change the METHOD_VERSION tripwire caught below); no floor/pair/divisor/clip/marker change.
# 58d3738 -> 8db37af (2026-07-24): upstream is now VENDORED as a submodule (external/Phenocycler_Analysis)
# instead of a working-tree checkout, so this SHA pin is authoritative rather than advisory. v7 -> v8:
# audited ADDITIVE / contract-safe -- adds phenocycler/reference_selection.py (the Gate-2 reference>=target
# frequency rule + ACCEPTED_REFERENCE_OVERRIDES) and a frozen ACCEPTED_PAIRS tuple with import-time
# validation. The ONLY new raises are the ACCEPTED_PAIRS-in-CANDIDATE_PAIRS membership check and a
# min_reference_to_target_ratio>0 guard. No change to MARKER_INPUT_FLOOR_METHODS, pair_input_qc,
# canonical_divisor, or input_qc_upper_quantile -- so CODEX's accepted first-pass divisors and 88/90
# screen states are unaffected. (The level-2 pair-web redesign will consume reference_selection.)
UPSTREAM_PINNED_SHA = "8db37af"

# METHOD_VERSION is the finer tripwire: the SHA pin above only catches COMMITTED moves, but upstream
# is developed with UNCOMMITTED working-tree edits, and METHOD_VERSION bumps ONLY on a real method
# change -- so this catches those even before they are committed. v6 -> v7 (2026-07-23, working tree
# ahead of 55b904d): the positive CALL is now divisor-gated only (`value > divisor`); the former
# NNMF-target-component rescue below the divisor is dropped (negligible for clean markers <=0.3 pp;
# it doubled high-background markers like CD11b). Divisors, floors, pairs and the clip are UNCHANGED,
# so CODEX's contract and its screen states/divisors are unaffected -- only per-cell calls + the
# review-overlay category labels change. When this trips, re-audit and bump.
# v7 -> v8 (2026-07-24, committed at 8db37af): adds Gate-2 reference selection + the frozen ACCEPTED_PAIRS
# set. Purely additive and contract-safe (audited above at UPSTREAM_PINNED_SHA); CODEX screen states and
# divisors are unaffected -- only new machinery the level-2 redesign opts into.
EXPECTED_METHOD_VERSION = "restore-pair-validation-v8"

# --- Level-2 sub-typing pairs (authoritative spec, user 2026-07-24) -----------------------------
# RESTORE thresholds are computed GLOBALLY per scene (whole-cell), ONE per target -- NOT within gated
# subsets. "Level 2" = additional targets appended to the pair list under the same counterpart rules;
# the compartment hierarchy is applied AFTERWARD (codex_lineage_run.py, Phase I) to INTERPRET the
# thresholds. This RETIRES the old within-compartment `IMMUNE_PAIRS` (CD3e<-CD68, CD4<-CD31, CD8a<-CD68,
# ...): every immune target now takes ECAD, matching the whole-scene convention and upstream v8's
# ACCEPTED_PAIRS direction. The old co-receptor hazard (CD4/CD8a<-CD3e) is thereby moot, but the guard
# in no_rare_counterparts() is kept as a defensive invariant.
#
# Counterpart rules (carried over from level 1): epithelial-family target -> CD31; hormone -> other
# hormone; neural -> panCK (epCAM dropped -- 8/15-donor only, and panCK measured a cleaner exocrine
# negative, see docs/codex_restore_pairs.md); LGALS3 -> CD31 (stated exception); everything else -> ECAD.
# Three markers do DOUBLE DUTY (PDPN, CD34, LGALS3): ONE threshold each, disambiguated by which parent
# node the cell landed in (the CD31 gate separates lymphatic-EC from fibroblast for PDPN/CD34). CD44,
# LGALS3 are STATE annotations, never a class call.
#
# 22 pairs. With FIRST_PASS_PAIRS(5) + HORMONE_COMPARATOR_PAIRS(3) the full screened set is 30 pairs
# over 27 markers (= PANEL_31 minus BActin/NaKATPase/BCatenin/ColIV; epCAM/KRT8-18 are 8/15-donor and
# also excluded). CHGA<-CD31 is shared with ENDOCRINE_SPLIT_PAIR (already reviewer-accepted, divisor ~650).
LEVEL2_PAIRS = [
    # === IMMUNE (CD45+) -> ECAD ===
    ["CD3e",   "ECAD"],   # T
    ["CD68",   "ECAD"],   # macrophage
    ["CD163",  "ECAD"],   # M2 polarization -- nested under CD68 at interpretation
    ["CD4",    "ECAD"],   # T-helper -- read ONLY within CD3e+ (monocytes/macrophages are CD4-dim)
    ["CD8a",   "ECAD"],   # T-cytotoxic
    ["HLADR",  "ECAD"],   # APC/activation -- WEAKEST pair (ductal upregulates HLA-II in T1D -> flag conservative)
    ["CD20",   "ECAD"],   # B -- re-included (CD20 unwind); prevalence-collapse marker
    ["CD35",   "ECAD"],   # FDC/TLS -- near-zero without lymphoid aggregates; prevalence-collapse marker
    # === EPITHELIAL (ECAD+): family -> CD31, endocrine hormones -> other hormone ===
    ["CHGA",   "CD31"],   # pan-endocrine split (== ENDOCRINE_SPLIT_PAIR; already VALID, divisor ~650)
    ["CK19",   "CD31"],   # ductal
    ["panCK",  "CD31"],   # epithelial confirmatory (also the neural counterpart)
    ["KRT14",  "CD31"],   # basal/metaplastic -- near-zero in normal tissue; prevalence-collapse marker
    ["INS",    "GCG"],    # beta
    ["GCG",    "INS"],    # alpha
    ["SST",    "INS"],    # delta -- sparse; prevalence-collapse marker
    # === ENDOTHELIAL (CD31+) -> ECAD (dual-use PDPN/CD34; CD31 gate disambiguates at interpretation) ===
    ["PDPN",   "ECAD"],   # NOT CD31 -- lymphatic EC are CD31+; PDPN+ -> lymphatic under Endothelial, CAF under Mesenchymal
    ["CD34",   "ECAD"],   # capillary EC / resting fibroblast (parent-resolved)
    # === NEURAL (PGP9.5+) -> panCK (epCAM substitute -- cohort-wide + cleaner negative) ===
    ["GAP43",  "panCK"],  # growth-associated/sprouting axons; prevalence-collapse marker
    ["B3TUBB", "panCK"],  # consensus pan-neural with PGP9.5
    # === MESENCHYMAL (VIM+) -> ECAD ===
    ["SMA",    "ECAD"],   # smooth muscle / pericyte / myofibroblast
    # === STATE (annotation, not a lineage class call) ===
    ["CD44",   "ECAD"],   # activation/stemness state
    ["LGALS3", "CD31"],   # dual macrophage+ductal -> ECAD disqualified; CD31 (stated exception). State only.
]

# Exocrine comparator references for the three hormones (user decision 2026-07-24). The hormone->hormone
# reference goes UNDERSIZED in exactly the disease states of interest (alpha<beta in ND -> INS<-GCG
# undersized; beta-loss in T1D -> GCG<-INS / SST<-INS undersized). These extra CANDIDATE pairs give the
# reference_selection dossier (Phase F) an exocrine comparator to weigh against the undersized hormone
# reference, instead of blindly accepting it via override. The hormone-hormone pair stays the predeclared
# baseline; these are additional candidates only (one-ref-per-target is required only for the frozen exit).
HORMONE_COMPARATOR_PAIRS = [
    ["INS", "panCK"],
    ["GCG", "panCK"],
    ["SST", "panCK"],
]

# All pairs the CODEX screen fits (adapters derive SCREENED_MARKERS / CODEX_PAIRS from this).
ALL_SCREEN_PAIRS = FIRST_PASS_PAIRS + LEVEL2_PAIRS + HORMONE_COMPARATOR_PAIRS

# Immune sub-split precedence, LOWEST first (later writes win, matching `lineage._cell_types`).
# B sits at the bottom deliberately: the abundant, well-supported types claim their cells before a
# CD20 call is allowed to stand, so a CD20+CD68+ cell resolves to Macrophage rather than B. This
# is already how upstream orders it (lineage.py:115-122) -- recorded here so a CODEX-specific
# rewrite of `_cell_types` cannot silently reintroduce the inverted order.
#
# B is now a SCREENED call (2026-07-24): CD20 is re-included as a RESTORE target (`CD20 <- ECAD`), so a
# `CD20_pos` column IS produced and the B gate can fire -- but only INSIDE the Immune (CD45+) compartment
# and at the BOTTOM of this precedence, so the abundant, well-supported types (Macrophage, T) claim their
# cells first and a CD20+CD68+ cell resolves to Macrophage. This bottom-precedence placement is what keeps
# CD20's known stray-on-epithelium signal (48.6% of CD20-only cells are ECAD+, see CD20_STRAY_SIGNAL) from
# producing spurious B calls: those cells are Epithelial, never reach the Immune sub-split. CODEX still
# lacks CD79a (Senior's B-cell fallback), so B rests on CD20 alone.
# DC / Neutrophil / Plasma / NK are unreachable in the 31-plex (no CD11c, CD209, MPO, CD38,
# CD56/CD57) and are omitted rather than left to evaluate as all-False.
IMMUNE_PRECEDENCE = ["Immune-other", "B", "Macrophage", "T"]
IMMUNE_SUBTYPE_GATES = {
    "B":          ["CD20"],   # screened (CD20 <- ECAD); called only within CD45+ at bottom precedence
    "Macrophage": ["CD68", "CD163"],
    "T":          ["CD3e"],
}

# CD20 IS NOT A COMPARTMENT ANCHOR, matching upstream 6d49ae7 which removed it from Senior's
# any-positive Immune union: a lone CD20 call is far likelier to be stray signal on an abundant
# cell than a real B cell, and a union would let it pull that cell into Immune before a more
# frequent compartment claimed it. CODEX is compliant by construction -- the Immune gate is CD45
# alone -- and B is reachable only as a sub-split INSIDE Immune, i.e. it already requires CD45+.
#
# Measured on the corrected CODEX data, and the margin is stark. Of 42,183 CD20+ cells at the
# donor-local floor, only 2,450 (5.8%) are also CD45+. 39,197 (92.9%) are CD20-only with no
# CD45, CD3e or CD68 -- and 48.6% of THOSE are ECAD+, i.e. abundant epithelium carrying stray
# CD20. So "CD20+ = 1.58% of cells" was measuring mostly noise; the CD45-gated figure of ~0.09%
# is a plausible pancreatic B-cell rate.
#
# Senior notes B cells still reach Immune via CD79a, which this panel lacks -- so with CD20 excluded
# from RESTORE (92f379a), B is an abstention (see IMMUNE_SUBTYPE_GATES above), not an Epithelial
# reroute. The measurement below is KEPT as the evidence that CD20 is untrustworthy at the single-cell
# level (92.9% of CD20+ cells are CD45-negative, half of those ECAD+ epithelium), which is exactly
# why upstream deferred it -- but it no longer supports "CD20 stays a target": there is no CD20 call.
CD20_NOT_A_COMPARTMENT_ANCHOR = True
CD20_STRAY_SIGNAL = {           # measured post-REDSEA at donor-local floors; retained as evidence only
    "cd20_pos": 42183,
    "cd20_and_cd45_pos": 2450,          # 5.8% -- the plausible B-cell population
    "cd20_only": 39197,                 # 92.9%, of which 48.6% are ECAD+
    "cd20_only_ecad_pos_pct": 48.6,
}

# The CD8a/HLADR "prevalence anomaly" -- RESOLVED 2026-07-23, and it was NOT a marker property.
#
# MECHANISM: upstream `raw_otsu_floor` CLIPS at the 99.5th percentile before running Otsu
# (restore_validation.py:772-804; the clip quantile is `input_qc_upper_quantile`, default 0.995). On
# the 2-3 donors where a sparse marker's positive rate is BELOW 0.5% (CD8a on 112/6356, HLADR on
# 6516/6356/6505, CHGA on 6548/6550/6551), the entire positive tail sits above p99.5, so the clip
# deletes it and Otsu bisects background. The bright cells are still there (6356 CD8a p99.9 = 23,792,
# max 48,684) -- the clipped threshold just cannot see them, landing at ~median instead of ~14,000.
#
# The distributions are NORMAL (cohort dynamic range p99.5/median: CD8a 70, HLADR 31, HIGHER than
# CD68's 35). The recorded 17.66/36.54 figures were the REDSEA-aggravated ones; the operative RAW
# figures were 6.6%/7.4%, concentrated on those 2-3 donors.
#
# FIX (validated on the actual RESTORE fit, not just the floor): screen the immune pairs UNCLIPPED
# (`input_qc_upper_quantile = 1.0`). Measured donor 6356: CD8a<-CD68 MODEL_UNSTABLE (floor 350 at
# median, divisor 349) -> VALID (floor 10,671 in the tail, divisor 1,343); HLADR<-CD31 MODEL_UNSTABLE
# -> VALID (floor 4,735). The good donor 6534 is unaffected (VALID both ways). Unclipped Otsu gives
# CD8a 0.6% / HLADR 0.5% cohort, sane and stable on all 15 donors, with CD8a+/CD3e+ ~= 0.95.
# `max(Otsu, Triangle)` was a REJECTED band-aid (a disguised p98-99 quantile cut). So the immune
# markers stay LINEAR_OTSU; the CLIP -- not the statistic -- is the control. The first-pass pairs keep
# the locked 0.995 clip (they were never broken; unclipping only perturbs their divisors).
#
# The two input-QC clips. The clip is a WHOLE-CONFIG scalar in upstream (`PairValidationConfig.
# input_qc_upper_quantile`, one float per evaluate_locked_cohort call), so a mixed-clip level-2 screen
# MUST run as TWO cohort passes (see codex_restore_validation): the 0.995 group for abundant/structural
# targets, and the 1.0 (unclipped) group for the clip-sensitive targets whose positive tail sits above
# p99.5 on low-rate donors.
LOCKED_UPPER_QUANTILE = 0.995     # the clipped group (structural + abundant + CHGA)
IMMUNE_SCREEN_UPPER_QUANTILE = 1.0  # the unclipped group value (kept name for back-compat)

# Targets screened UNCLIPPED. Membership = the codex_floors.clip_sensitive set MINUS CHGA (CHGA rides
# the 0.995 group where its max()@0.995 divisor is robust). A pair goes in the unclipped group iff its
# TARGET is here; otherwise the 0.995 group. MEASURED by `codex_floors.py` (2026-07-24, all 15 donors
# RAW) -- this literal is the single source of truth for the clip split.
#
# These 10 markers have min(p99.5/median) < 7 AND prev_max < 10%: their bright positive population sits
# in the tail ABOVE the p99.5 clip on the low-rate donors, so a clipped 2-component Otsu bisects
# background. Note two SPARSE markers are NOT here: CD20 (prev 14.5% at its low stray-signal floor ->
# abundant-at-floor, not clipped-out) and GAP43 (min_dynrange 13.6 -> long tail, clip is safe). Both
# still ride the borrow/global-divisor prevalence-collapse policy (PREVALENCE_COLLAPSE_MARKERS), which
# is a positivity-step concern, not a clip concern.
#
# PDPN was REMOVED from this set (2026-07-24, image review): the clip_sensitive HEURISTIC flagged it
# (min_dynrange 2.2), but a clip sweep proved unclipping is the WORST choice for it -- unclipped gives
# 0.18% median called with 3 zero-donors + degeneracy, while 0.995 gives 0.43% median (up to 5.2%),
# all 15 donors stable. PDPN's positives are NOT hidden above the clip (the CD8a/HLADR mechanism);
# unclipping just raises its divisor and starves the fit. So PDPN rides the 0.995 group. (Two donors
# read ~0 PDPN at every clip -- genuinely PDPN-absent scenes, handled by the prevalence policy.)
CLIP_UNCLIPPED_TARGETS = frozenset({
    "CD3e", "CD4", "CD8a", "HLADR",           # lymphocyte/activation
    "CD35", "KRT14", "SST", "INS",            # sparse subtype targets (INS bright, sparse cohort-wide)
    "CD44", "LGALS3",                         # low-dynamic-range state/dual-use targets
})

# Targets whose POSITIVE mode collapses on some donors: RESTORE's 2-component fit then splits the
# NEGATIVE distribution instead of failing, and the `counterpart > median` prefilter selects half the
# tissue while the target positive mode may be < 0.1%. This is DISTINCT from reference sparsity
# (reference_selection's job) and from the p99.5 clip artifact (CLIP_UNCLIPPED_TARGETS). Policy
# (user 2026-07-24): estimate the threshold on scenes KNOWN to contain positives and apply it GLOBALLY,
# and/or check positive-cluster size before trusting any per-scene value. Consumed by
# codex_positivity.resolve_divisors (Phase H): a per-donor divisor is trusted only if its positive
# cluster is adequate, else a cohort divisor (median of adequate donors) is borrowed; abstain
# cohort-wide only if NO donor has an adequate positive cluster.
PREVALENCE_COLLAPSE_MARKERS = ("KRT14", "CD35", "CD20", "SST", "GAP43")

CLIP_ARTIFACT_DONORS = {          # marker: donors where the p99.5 clip broke plain Otsu (RAW)
    "CD8a":  ["112", "6356"],
    "HLADR": ["6516", "6356", "6505"],
    "CHGA":  ["6548", "6550", "6551"],
}

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
# docs/codex_restore_pairs.md). Triangle falls below Otsu for all five STRUCTURAL gate markers
# (CD45/ECAD/CD31/PGP9.5/VIM), so max(Otsu, Triangle) reduces to plain Otsu there and the second
# statistic buys nothing. CHGA is the exception -- Triangle EXCEEDS Otsu on 9/15 donors (RAW),
# tightening cohort-mean 5.30% -> 1.50% high. This "Triangle < Otsu" claim holds ONLY for the six
# abundant structural markers; it does NOT generalise to the sparse immune markers (Triangle > Otsu
# for CD8a 2/15, CD3e 6/15, HLADR 7/15 on RAW), which is exactly why upstream keeps CD3e/CD20 on
# max() and why floor_method() fails closed rather than letting Senior's map be inherited blindly.
LINEAR_OTSU_FLOOR = "linear_otsu"
MAX_LINEAR_OTSU_TRIANGLE_FLOOR = "max_linear_otsu_triangle"

# Measured by `codex_floors.py` on the RAW parquet, all 15 donors. The rule (codex_floors.choose_method):
# max(Otsu,Triangle) iff a marker is BOTH low-dynamic-range (min p99.5/median < 7) AND sparse
# (< 10% at the max floor); plain Otsu otherwise. `min_dynrange` is the per-donor minimum of
# p99.5/median -- the predictor of the p99.5-clip failure. Re-derive with `python codex_floors.py`.
MARKER_INPUT_FLOOR_METHODS = {
    # --- structural gates: high dynamic range, plain Otsu is safe -----------------------------
    "CD45":   LINEAR_OTSU_FLOOR,             # min_dynrange 32.5, prev 2.97%
    "ECAD":   LINEAR_OTSU_FLOOR,             # min_dynrange  4.0 but ABUNDANT (31.25%) -> Otsu, not max
    "CD31":   LINEAR_OTSU_FLOOR,             # min_dynrange 55.7, prev 4.56%
    "PGP9.5": LINEAR_OTSU_FLOOR,             # min_dynrange 83.3, prev 2.05%
    "VIM":    LINEAR_OTSU_FLOOR,             # min_dynrange 97.7, prev 4.89%
    "panCK":  LINEAR_OTSU_FLOOR,             # min_dynrange 15.5, ABUNDANT (15.73%)
    # --- myeloid: sparse but high dynamic range, plain Otsu is safe at 0.995 -------------------
    "CD68":   LINEAR_OTSU_FLOOR,             # min_dynrange 13.2, prev 2.09%
    "CD163":  LINEAR_OTSU_FLOOR,             # min_dynrange 22.7, prev 1.79%
    # --- lymphocyte/activation: plain Otsu, but screened UNCLIPPED (IMMUNE_SCREEN_UPPER_QUANTILE).
    #     At the locked 0.995 clip these break on the 2-3 low-rate donors; unclipped they are clean on
    #     all 15 (CD8a 0.6%, HLADR 0.5%, CD3e 0.6%, CD4 0.8%) and the RESTORE fit goes VALID. max() was
    #     rejected -- see the CD8a/HLADR resolution note above. CHGA stays max() because it rides the
    #     first-pass run at the locked 0.995 clip, where its divisor is robust.
    "CHGA":   MAX_LINEAR_OTSU_TRIANGLE_FLOOR,  # first-pass @ 0.995; Triangle wins 9/15, divisor robust
    "CD3e":   LINEAR_OTSU_FLOOR,             # immune @ unclipped -> 0.6% cohort, sane on all 15
    "CD4":    LINEAR_OTSU_FLOOR,             # immune @ unclipped -> 0.8%
    "CD8a":   LINEAR_OTSU_FLOOR,             # immune @ unclipped -> 0.6%, CD8a+/CD3e+ ~= 0.95
    "HLADR":  LINEAR_OTSU_FLOOR,             # immune @ unclipped -> 0.5%
    # --- level-2 sub-typing markers (2026-07-24). ALL plain Otsu: max(Otsu,Triangle) was rejected
    #     (a disguised quantile cut); the sparse ones ride the UNCLIPPED group (CLIP_UNCLIPPED_TARGETS)
    #     instead of a max() band-aid. Floors measured on RAW by codex_floors.py; see OBSERVED_FLOORS_RAW.
    "CD20":   LINEAR_OTSU_FLOOR,             # B -- sparse, unclipped; prevalence-collapse
    "CD35":   LINEAR_OTSU_FLOOR,             # FDC/TLS -- sparse, unclipped; prevalence-collapse
    "SST":    LINEAR_OTSU_FLOOR,             # delta hormone -- sparse, unclipped; prevalence-collapse
    "GAP43":  LINEAR_OTSU_FLOOR,             # neural -- sparse, unclipped; prevalence-collapse
    "KRT14":  LINEAR_OTSU_FLOOR,             # basal/metaplastic -- sparse, unclipped; prevalence-collapse
    "INS":    LINEAR_OTSU_FLOOR,             # beta hormone (abundant in ND, sparse in T1D)
    "GCG":    LINEAR_OTSU_FLOOR,             # alpha hormone
    "CK19":   LINEAR_OTSU_FLOOR,             # ductal epithelium
    "B3TUBB": LINEAR_OTSU_FLOOR,             # pan-neural consensus
    "CD34":   LINEAR_OTSU_FLOOR,             # capillary EC / resting fibroblast (dual-use)
    "PDPN":   LINEAR_OTSU_FLOOR,             # lymphatic EC / inflammatory-CAF (dual-use)
    "SMA":    LINEAR_OTSU_FLOOR,             # smooth muscle / pericyte
    "CD44":   LINEAR_OTSU_FLOOR,             # activation/stemness state
    "LGALS3": LINEAR_OTSU_FLOOR,             # macrophage+ductal state
}

# CORRECTED 2026-07-23: RESTORE floors on the RAW parquet, NOT the REDSEA one. Upstream's
# `pair_input_qc(target_raw, reference_raw, ...)` (restore_validation.py:1906) floors on
# `raw__{marker}`, and `load_validation_sample` (:1499-1512) np.allclose-asserts that `raw__{m}`
# equals the UNCORRECTED QuPath `Cell: {m}: Mean`. REDSEA enters only as the divisor/separator axis
# (`redsea__{m}`). So `OBSERVED_FLOORS_RAW` below is the OPERATIVE table; `OBSERVED_FLOORS`
# (post-REDSEA) is reference-only. The earlier note here had this exactly backwards, and the
# HANDOFF's "measure floors post-REDSEA" instruction was wrong on the same premise.
#
# The post-REDSEA table is still recorded (correction lowers every floor, -3.8% CD31 to -40% CHGA)
# as a diagnostic, but a rerun should be regression-checked against OBSERVED_FLOORS_RAW.
OBSERVED_FLOORS = {          # marker: (median, pct_high_at_median, spread_x) -- post-REDSEA, REFERENCE ONLY
    "CD45":   (2641.3,  1.91,  3.6),
    "ECAD":   (9213.8, 29.24,  2.0),
    "CD31":   (2968.0,  4.16,  4.1),
    "PGP9.5": (4658.8,  1.62, 10.1),
    "VIM":    (4176.3,  4.04,  2.5),
    "panCK":  (3308.8, 14.05,  2.6),
    "CHGA":    (811.9,  1.33,  4.3),   # was (570.8,1.23,11.4): a chimera taking floor+spread from
                                       # plain Otsu while CHGA is registered as max(Otsu,Triangle).
                                       # Corrected to the declared statistic (measured 2026-07-23).
}
OBSERVED_FLOORS_RAW = {      # pre-REDSEA -- THIS is what RESTORE floors on. Max-statistic for CHGA only.
    # Measured 2026-07-24 by codex_floors.py on the RAW parquet (median across 15 donors, each marker's
    # registered statistic). All 27 screened markers are recorded so `--check` regresses the full set.
    # --- level-1 structural gates + CHGA ---
    "CD45": 2871.7, "ECAD": 12622.9, "CD31": 3084.9, "PGP9.5": 5547.0,
    "VIM": 4459.5, "panCK": 4092.7, "CHGA": 950.9,   # CHGA = max(Otsu,Triangle); panCK median of 15
    # --- immune sub-typing ---
    "CD3e": 540.5, "CD4": 746.6, "CD8a": 5256.5, "CD68": 716.8, "CD163": 1240.1,
    "HLADR": 818.1, "CD20": 279.1, "CD35": 184.9,
    # --- epithelial family + hormones ---
    "CK19": 6719.1, "KRT14": 2069.9, "INS": 203.0, "GCG": 3617.9, "SST": 201.7,
    # --- endothelial (dual-use) ---
    "PDPN": 712.5, "CD34": 2603.3,
    # --- neural ---
    "GAP43": 855.7, "B3TUBB": 2273.9,
    # --- mesenchymal + state ---
    "SMA": 7381.1, "CD44": 1104.5, "LGALS3": 10065.5,
}

# Fraction of cells REDSEA clamped to ZERO, cohort median (max donor in parentheses).
#
# RESOLVED 2026-07-23 (handoff issue 3, CLOSED, not a blocker): the clamp lands on BACKGROUND, not
# on gate-positive cells. Mesenchymal moves only 2.95% -> 2.79% (-0.16 pp); only ~1.1% of raw-VIM+
# cells clamp to zero, and the VIM-positive set has Jaccard 84-88% raw-vs-corrected. Decisive
# counterfactual: with prevalence-matched corrected floors, pre/post Other differs by 0.15 pp and
# Epithelial by 0.00 pp -- REDSEA does not change WHICH cells rank positive; the apparent +4.5 pp
# Other shift is just the Otsu cut re-landing at a higher percentile of the corrected distribution.
# It is moot for RESTORE anyway, which floors on RAW (see OBSERVED_FLOORS_RAW). Kept as a diagnostic.
REDSEA_ZERO_CLAMP = {        # marker: (cohort_median_pct, worst_pct, worst_donor)
    "CD45":   (16.6, 32.3, "6534"),
    "ECAD":   (31.1, 35.8, "6505"),
    "CD31":   (20.8, 28.4, "6450"),
    "PGP9.5": (17.4, 29.7, "6450"),
    "VIM":    (39.7, 54.6, "6534"),
    "panCK":  (28.9, 39.6, "6534"),
    "CHGA":    (4.5,  6.8, "6549"),
}

# A hard ceiling discovered while closing issue 3: post-REDSEA only 69.03% of cells have ECAD > 0 at
# all (per-donor 64.2%-80.5%), versus 99.97% on the raw parquet -- subtract-only over-subtraction on
# the abundant acinar mass. No positive ECAD threshold can therefore assign more than ~69% of cells
# to Epithelial on cells_redsea, and the five-gate nonzero union caps non-Other at 98.60%. This is
# one more reason RESTORE (and its floors) run on RAW, where ECAD > 0 on essentially every cell.
ECAD_NONZERO_POST_REDSEA = 0.6903
FIVE_GATE_NONZERO_UNION_POST_REDSEA = 0.9860


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
        "CORRECTED 2026-07-23: at the donor-local FLOORS (a proxy, NOT RESTORE positivity), Other is "
        "the LARGEST compartment, ~61% post-REDSEA / ~56% raw, and ~62% of it is app-labelled Acinar "
        "(61% of all Acinar fails every gate). Epithelial lands at ~28% in a tissue the app puts at "
        "~72% epithelial, and 67% of app-Beta cells route to Neural. This is expected of Otsu FLOORS "
        "used as positivity calls -- the post-REDSEA ECAD floor (9214) sits ABOVE the raw ECAD median, "
        "so Otsu splits the positive mode. It is NOT the acceptance verdict: RESTORE's per-image "
        "divisor is what supplies positivity (and is measurably more permissive -- donor 6516 ECAD "
        "divisor 9023 vs input floor 11834). The composition MUST be re-judged on RESTORE output "
        "(plan Phase 5); if Other stays ~60% there, the gate marker set is reopened as measured "
        "evidence. The old 'Other is small (adipocytes/Schwann/stellate)' claim was never measured."
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


def no_rare_counterparts() -> None:
    """No pair may use an EXCLUDED marker in either role, or be a self-pair.

    Assertions rather than comments because a bad reference yields plausible output rather than
    an error. Mirrors upstream's both-role import guard in `restore_validation` (92f379a); kept
    here because the CODEX pair sets never pass through `CANDIDATE_PAIRS`.
    """
    for tgt, cp in ALL_SCREEN_PAIRS:
        for m, role in ((tgt, "target"), (cp, "reference")):
            if m in RESTORE_EXCLUDED_MARKERS:
                raise ValueError(
                    f"{tgt} <- {cp}: {m} is excluded from RESTORE in both roles (here as the "
                    f"{role}). {RESTORE_EXCLUDED_MARKERS[m]}"
                )
        if tgt == cp:
            raise ValueError(f"{tgt} cannot be its own reference")
    # CD4/CD8a are T-cell co-receptors: CD3e is definitionally NOT exclusive of them. Vacuous now that
    # every immune target takes ECAD, but kept as a defensive invariant against a future re-edit.
    for tgt, cp in LEVEL2_PAIRS:
        if tgt in ("CD4", "CD8a") and cp == "CD3e":
            raise ValueError(
                f"{tgt} <- CD3e is biologically invalid: {tgt} is a T-cell co-receptor and is "
                "co-expressed with CD3e. Use a myeloid (CD68/CD163) or non-immune (CD31) counterpart."
            )


def assert_pairs_valid(marker_pairs) -> None:
    """Validate a pair list at the point of use, not just at import.

    UPSTREAM CONTEXT (2026-07-23): `config.MARKER_PAIRS` / `FUNCTIONAL_PAIRS` and the
    `run_restore` / `run_restore_functional` entry points were RETIRED rather than repaired --
    seven MARKER_PAIRS and three FUNCTIONAL_PAIRS referenced CD20, and upstream's note states
    there is no mechanical substitute because "CD3e is invalid for the T-subset targets
    (CD4/CD8/FOXP3 cells are CD3e+)". That is the same objection measured here (CD8a<-CD3e
    58% overlap by the p99 proxy; co-expression by biology).

    So there is no longer an implicit CD20-referenced default to fall back into -- the hazard
    that motivated the earlier guard is gone at the source. This check remains because CODEX
    pair sets never pass through upstream's `CANDIDATE_PAIRS` import guard, so they would
    otherwise be validated nowhere.
    """
    if not marker_pairs:
        raise ValueError("no marker_pairs supplied; CODEX pairs must always be passed explicitly")
    for tgt, cp in marker_pairs:
        for m, role in ((tgt, "target"), (cp, "reference")):
            if m in RESTORE_EXCLUDED_MARKERS:
                raise ValueError(
                    f"{tgt} <- {cp}: {m} is excluded from RESTORE in both roles (here as the "
                    f"{role}). {RESTORE_EXCLUDED_MARKERS[m]}"
                )
        if tgt == cp:
            raise ValueError(f"{tgt} cannot be its own reference")


def validate_against_panel() -> None:
    """Every marker referenced here must exist in the 31-marker panel, on every donor."""
    used = {SEGMENTATION_MARKER}
    for ms in COMPARTMENT_GATES.values():
        used.update(ms)
    for ms in IMMUNE_SUBTYPE_GATES.values():
        used.update(ms)
    used.update(NEURAL_CONSENSUS)
    used.update(MESENCHYMAL_SUBTYPE)
    for tgt, cp in ALL_SCREEN_PAIRS:
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
no_rare_counterparts()
