"""codex_typing.py -- the CODEX side of the Islet-Explorer-Senior cell typing port (owner, 2026-10-16).

The owner's direction: "use the cell typing from Islet-Explorer-Senior as much as possible" (D3, D4, D5), and
"drop cell classes that indicate neural, blood vessel or lymphatic and use the objects instead" (A6).

Senior's typing is the v3 worktree `~/IO60panc2nd/Islet-Explorer-Senior-v3`, the run of 2026-10-08
(`/data/io60-calibration/own_compartment_frozen_reference_repair_20261008/`), whose source snapshots are
byte-identical to the worktree files listed in `../PORTED_FROM.json` under "typing". Its chain:

  REDSEA spillover correction, one compartment per marker  -> redsea_cells.py
  UniFORM rigid log-shift per (donor, marker); cohort floor/top scale to [0, 1]  -> align_scale.py
  context: islet object, structure object (Blood vessel > Nerve > Lymphatic) or ordinary tissue  -> cell_context.py
  structure, endocrine, immune-definition, supporting-tier, Exocrine-residual calls  -> classify.py
  refinement (neighbour concordance, own-marker bound), applied once  -> refine.py
  immune class and subtype, most specific rule wins  -> immune_subtype.py

Everything below is a fact about THIS panel that the Senior code reads from its own registry; every departure
from Senior is commented where it is made and repeated in README.md.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNIFORM = HERE.parent
REPO = UNIFORM.parents[1]
COHORT = UNIFORM / 'cohort.csv'
PIXEL_UM = 0.3774

# Senior marker name -> this panel's channel name, for every Senior marker this panel carries. Senior's
# M2Gal3 is the Mac-2 antibody, which recognises galectin-3, the protein the CODEX LGALS3 channel stains.
SENIOR_TO_CODEX: dict[str, str] = {
    'INS': 'INS', 'GCG': 'GCG', 'SST': 'SST', 'B3TUBB': 'B3TUBB', 'Podoplanin': 'PDPN', 'CD34': 'CD34',
    'CD31': 'CD31', 'E_cadherin': 'ECAD', 'Pan_Cytokeratin': 'panCK', 'CD3e': 'CD3e', 'CD4': 'CD4',
    'CD8': 'CD8a', 'CD20': 'CD20', 'CD68': 'CD68', 'CD163': 'CD163', 'HLA_DR': 'HLADR', 'SMA': 'SMA',
    'Vimentin': 'VIM', 'Collagen_IV': 'ColIV', 'CD44': 'CD44', 'Beta_actin': 'BActin',
    'b_Catenin1': 'BCatenin', 'Ker8_18': 'KRT8-18', 'EpCAM': 'epCAM', 'M2Gal3': 'LGALS3', 'Ki67': 'Ki67',
    # CD99 has no counterpart; CHGA takes its role (owner, 2026-10-15), as in the islet fit.
    'CD99': 'CHGA',
}

# Senior markers the typing reads that this panel never acquired. A rule term on one of them is removed and
# the removal recorded (Senior's own wording, celltype/score.meets_requirements: "a caller that means to
# waive a constraint on an unmeasured marker removes it from reqs and records that it did").
NOT_ON_PANEL = ('CD56', 'Caveolin', 'CD79a', 'CD57', 'Granzyme_B', 'CD11b', 'CD11c', 'CD209', 'CD66',
                'Iba1', 'MPO', 'IDO1', 'CD38', 'FOXP3', 'ICOS', 'CD39', 'PD_1', 'LAG3', 'TOX', 'VISTA',
                'TCF_1', 'Bcl_2', 'CD107a', 'iNOS', 'PD_L1', 'CD206', 'Keratin_5', 'TP63', 'SOX2', 'PCNA',
                'HLA_A')

# Acquired on 13 of the 20 images only (Aab_6521, Aab_6549, ND_0112, ND_6548, T1D_6533, T1D_6534 and
# T1D_6551 lack both). Senior's rule test reads an unmeasured negative as UNKNOWN, so every immune definition
# (all of them forbid EpCAM and Ker8_18) could never fire on those seven donors. They are left out of the
# typing on every donor, so all 20 are asked the same question.
NOT_ON_EVERY_IMAGE = ('KRT8-18', 'epCAM', 'Ghrelin')

# The compartment each channel is measured in: Senior's ruling
# (results/v3/release_20261003/calibration_review/marker_redsea_ruling.csv) for the markers it rules on, and
# the same reasoning for this panel's own markers (secretory granule and cytoskeletal markers Cytoplasm, like
# the hormones and keratins; surface markers Membrane, like every CD marker).
COMPARTMENT: dict[str, str] = {
    # Senior's ruling
    'INS': 'Cytoplasm', 'GCG': 'Cytoplasm', 'SST': 'Cytoplasm', 'B3TUBB': 'Cytoplasm', 'panCK': 'Cytoplasm',
    'VIM': 'Cytoplasm', 'SMA': 'Cytoplasm', 'ColIV': 'Cytoplasm', 'BActin': 'Cytoplasm', 'LGALS3': 'Cytoplasm',
    'KRT8-18': 'Cytoplasm', 'CD31': 'Membrane', 'CD8a': 'Membrane', 'CD3e': 'Membrane', 'CD4': 'Membrane',
    'CD34': 'Membrane', 'HLADR': 'Membrane', 'PDPN': 'Membrane', 'ECAD': 'Membrane', 'CD163': 'Membrane',
    'CD20': 'Membrane', 'CD44': 'Membrane', 'CD68': 'Membrane', 'BCatenin': 'Membrane', 'epCAM': 'Membrane',
    'Ki67': 'Nucleus',
    # this panel's own markers, by analogy
    'CHGA': 'Cytoplasm', 'PGP9.5': 'Cytoplasm', 'GAP43': 'Cytoplasm', 'CK19': 'Cytoplasm',
    'KRT14': 'Cytoplasm', 'Ghrelin': 'Cytoplasm', 'CD45': 'Membrane', 'CD35': 'Membrane',
    'NaKATPase': 'Membrane',
}
# Senior's marker registry: channels REDSEA passes through uncorrected (spillover_policy != redsea_corrected).
NO_SPILLOVER = ('DAPI', 'ColIV', 'Ki67')
# E-cadherin is measured twice (owner's acceptance criterion, 2026-10-26: a measurement is judged by its effect on
# positivity; measured 2026-10-28, output/positivity_method_comparison_20261028/redsea_zeroing_by_marker.csv).
# E-cadherin is a junction shared by adjacent epithelial cells, so REDSEA's spillover compensation subtracts a cell's
# own junction: it sets E-cadherin to exactly 0 in 28 % of cells, 6.4 % of all cells (13 % in ND_6232 and Aab_6090)
# being at the epithelial level uncompensated, and those acinar cells fell to the top-marker step (no other rule
# marker loses more than 0.3 % of cells to zeroing). Read uncompensated everywhere, the same junction becomes the
# neighbor's signal in an immune cell pressed against acinar cells: 17,884 of v16's rule-typed immune cells (9,479
# macrophages, 5,927 CD45+, 1,707 T, 730 B) failed the epithelial condition. So:
#   'ECAD'               REDSEA-compensated (Senior's ruling): the immune rules' epithelial condition, where the
#                        question is whether epithelial signal belongs to the cell rather than its neighbor;
#   'ECAD_uncompensated' QuPath's own membrane mean (`$OUT/cells/<image>.parquet`): the Exocrine rule and the top-marker
#                        step, where the question is whether an unclaimed cell has an E-cadherin junction at all.
UNCOMPENSATED_COPIES = {'ECAD': 'ECAD_uncompensated'}


def channel_of(marker: str) -> str:
    """The image channel a typing marker is measured on ('ECAD_uncompensated' -> 'ECAD')."""
    return next((m for m, c in UNCOMPENSATED_COPIES.items() if c == marker), marker)
NOT_FOR_TYPING = ('DAPI',)

# REDSEA parameters of the 2026-10-08 run (phenocycler/config.py).
REDSEA = {'comp_mode': 1, 'alpha': 1.0, 'edge_radius': 0, 'gap_bridge': 1, 'norm_form': 'donor'}

# Senior geometry QC, as its pipeline config sets it (phenocycler/config.py cell_qc_*): an InstanSeg object
# without a nucleus is excluded from analysis and from spillover context.
GEOMETRY_QC = {'min_cell_area_um2': 20.0, 'max_cell_area_median_multiple': 4.0, 'min_cell_solidity': 0.70,
               'min_raster_area_ratio': 0.50, 'max_raster_area_ratio': 2.00, 'no_cytoplasm_nc_ratio': 0.999,
               'duplicate_radius_um': 1.0, 'duplicate_area_tolerance': 0.15, 'missing_geometry_policy': 'exclude'}
# Which of Senior's QC flags EXCLUDE a cell here (owner, 2026-10-16: "Cell quality check should drop cells without
# a nucleus, and that should take care of small cells, but do not drop irregular cells"). Size, solidity and the
# raster/reported area ratio are measured and written in `qc_flags`, but they exclude nothing; what excludes is a
# missing nucleus and the technical failures that leave a cell with no measurement at all.
QC_EXCLUDES = ('missing_nucleus_geometry', 'nucleus_raster_empty', 'raster_empty', 'invalid_centroid',
               'invalid_cell_area', 'duplicate_detection')

# --- the classification (typing/subpanel.py, subpanel_io.py), in this panel's names --------------------------
HORMONES = ('INS', 'GCG', 'SST')
HORMONE_NAMES = {'INS': 'Beta', 'GCG': 'Alpha', 'SST': 'Delta'}
LEVEL_MARKER = 'CHGA'          # Senior's CD99
# Owner, 2026-10-18: "These hormones are mutually exclusive so there needs to be resolution or exclusion based on the
# approved islet cell typing methods ... We need to add a CHGA+INS-GCG-SST- cell type."
# Resolution: the strongest scaled hormone names the cell only when it is at least HORMONE_MARGIN times the
# runner-up (Senior's confirmation margin, islet_identity.REFERENCE_MARGIN). Below it the hormones cannot be told
# apart -- on this segmentation the runner-up is mostly a touching neighbour's hormone (it is 4-7x larger when that
# neighbour type surrounds the cell) -- and the cell is Excluded (HORMONE_UNRESOLVED_REASON).
HORMONE_MARGIN = 2.0
HORMONE_UNRESOLVED_REASON = 'hormone_unresolved_margin'
# Senior's Endo-hormone-negative (celltype/run_islet.py, run_sameslide_protein.py: inside an islet, hormone, then
# CD99, then everything else): inside an Islet or Cluster, CHGA at or above the image's ordinary-tissue CHGA level
# and every hormone measured and below the image's hormone level.
HORMONE_NEGATIVE = 'CHGA+ hormone-negative'
# Senior's structure classes, strongest first. The objects here are the original project's pixel-classifier
# annotations (Capillary on CD31, Nerve on B3TUBB and INS, Lymphatics on PDPN); Capillary is Senior's
# "Blood vessel".
STRUCTURES = ('Blood vessel', 'Nerve', 'Lymphatic')
OBJECT_CLASS = {'Blood vessel': 'Capillary', 'Nerve': 'Nerve', 'Lymphatic': 'Lymphatic'}
# The markers that answer each context's question (subpanel_log.json "subpanels"). Senior's rule is that a
# vessel is read on the marker that DREW its objects (Caveolin/CD34 there, "NOT CD31: ... CD31 was never what
# made them"); here CD31 is what drew them. CD56 and Caveolin are not on this panel.
SUBPANEL = {'Nerve': ('B3TUBB',), 'Lymphatic': ('PDPN',), 'Blood vessel': ('CD31',),
            'exocrine': ('ECAD', 'panCK')}
CLASSIFIER_MARKER = {'Nerve': ('B3TUBB',), 'Lymphatic': ('PDPN',), 'Blood vessel': ('CD31',)}
# Senior's bottom_tier_immune on this panel (CD206, Granzyme_B, CD107a, PD_1, TOX, TCF_1, ICOS, CD39 and
# iNOS are not on it).
BOTTOM_TIER_IMMUNE = ('CD4', 'CD8a', 'CD163', 'HLADR', 'LGALS3')
# Senior's lineages.LINEAGE_MARKERS on this panel: the identity markers whose highest expressed value decides
# correction C (a structure call needs the object's own marker on top). Keratin_5, TP63, Caveolin, CD56,
# Iba1, CD206, CD11c, CD209, CD11b, MPO, CD66, FOXP3 and CD79a are not on the panel; Ker8_18 and EpCAM are
# left out with NOT_ON_EVERY_IMAGE.
LINEAGE_MARKERS = {
    'Endocrine': ('INS', 'GCG', 'SST'), 'Exocrine': ('panCK', 'ECAD'),
    'Vascular': ('CD31', 'PDPN', 'CD34'), 'Stromal': ('SMA', 'VIM', 'ColIV'),
    'T_NK': ('CD3e', 'CD8a', 'CD4'), 'B_Plasma': ('CD20',), 'Myeloid': ('CD68', 'CD163', 'LGALS3'),
    'Neural': ('B3TUBB',),
}
SUPPORTING_ONLY = ('LGALS3',)   # Senior: M2Gal3, CD11b, FOXP3, CD56, TP63
IDENTITY_MARKERS = tuple(sorted({m for ms in LINEAGE_MARKERS.values() for m in ms}))
T_TIE_PROFILE_MARKERS = IDENTITY_MARKERS + ('CD45',)
# Owner, 2026-10-27: "The established cell typing should classify cells appropriately and there should be no default
# of any class such as exocrine or immune. If typing has 'left over' cells and a check of top markers can't place
# them, call them by their top marker e.g. CD45+". Exocrine is a rule (ECAD or panCK at IMMUNE_REQUIRE_MIN), Senior's
# supporting tier is gone, and a cell no rule claims is placed by its highest marker among TOP_MARKER_PROFILE
# (scaled, before the clip) when that marker is a class's defining marker, otherwise labeled '<marker>+'.
# the top-marker step asks whether an unclaimed cell has the marker at all: E-cadherin uncompensated (UNCOMPENSATED_COPIES)
TOP_MARKER_PROFILE = tuple(UNCOMPENSATED_COPIES.get(m, m) for m in T_TIE_PROFILE_MARKERS)
# Owner, 2026-10-27 ("Floor at zero"): for markers most cells express, the most common positive value is the expressing
# population, not background (REDSEA already compensates their background to 0, 17-28 % of cells); align_scale.py's
# floor sat above 71-86 % of cells (ECAD 34,035 against a median of 9,044). Their cohort-scale floor is 0.
ZERO_FLOOR_MARKERS = ('ECAD', 'KRT8-18', 'BCatenin', 'BActin', 'epCAM')
TOP_MARKER_PLACEMENT = {'CD3e': 'Imm-T-cell', 'CD20': 'Imm-B-cell', 'CD68': 'Imm-macrophage', 'CD163': 'Imm-macrophage',
                        'ECAD_uncompensated': 'Exocrine', 'panCK': 'Exocrine'}
# The top marker is taken at any level (owner: "Any level"). ColIV is extracellular matrix: "ColIV cannot be a cell
# type, if it is highest, it is a reason for exclusion".
TOP_MARKER_EXCLUDE = ('ColIV',)
TOP_MARKER_EXCLUDED_REASON = 'top_marker_ColIV'
# Owner, 2026-10-27: the CD45 rule's class is "CD45+" (not "Immune"); a left-over cell whose top marker is CD45 is kept
# apart from it ("Keep them separate") as "CD45+ (top marker)".
CD45_RULE_LABEL = 'CD45+'


def top_marker_label(m: str) -> str:
    base = channel_of(m)                         # 'ECAD_uncompensated' is named by its channel
    return f'{base}+ (top marker)' if f'{base}+' == CD45_RULE_LABEL else f'{base}+'


# The top-marker step uses the Xenium pipeline's canary-marker logic (owner, 2026-10-28: "The top marker needs to use
# the canary marker logic that the Xenium processing uses for canary genes"; Islet-Explorer-Senior
# analysis/xenium/notebooks/02_phenotyping.ipynb cell 22, scripts/tissue_markers/pancreas.py LINEAGE_CANARY_PANELS /
# LINEAGE_EXCLUSIVE_GROUP). For a cell no rule claims:
#   depth     per lineage panel, the number of its markers positive (above the image's own level; Xenium: distinct
#             panel genes at lognorm >= log1p(3)); sum = the panel markers' excess over their levels (Xenium: the
#             summed lognorm); lineages rank by depth x 1000 + sum
#   groups    lineages map to mutually exclusive groups; the number of groups with a detected lineage decides:
#             1 group: that group's top lineage; 2 groups: the top lineage when its depth is >= CANARY_MIN_DEPTH and
#             >= CANARY_MARGIN more than the runner-up's, else 'lineage_ambiguous'; >= 3 groups: 'doublet_suspected'
#             (Xenium tags both and keeps the cell's call; here both are kept, named by their top detected marker and
#             never placed in a class, with the status in resolved_by)
#   naming    the lineage names an approved class (CANARY_PLACEMENT) when its depth is >= CANARY_MIN_DEPTH (Xenium's
#             two-distinct-marker evidence; its thin-panel lesson) and the class's defining marker is among the
#             detected ones; otherwise the cell is called by its top detected marker in that lineage ('<marker>+')
#   none      a cell with no lineage detected is called by its top marker at any level (owner, 2026-10-27: "Any
#             level"; the marker furthest above its own level), never placed in a class
#   ColIV     a cell called by ColIV is Excluded (owner, 2026-10-27)
# The panels are Senior's LINEAGE_MARKERS on this panel, E-cadherin read uncompensated, with CD45 in the lymphoid
# panels (Xenium's T panel carries PTPRC; the approved B-cell rule takes CD45 for CD79a).
CANARY_PANELS = {lin: tuple(UNCOMPENSATED_COPIES.get(m, m) for m in ms) for lin, ms in LINEAGE_MARKERS.items()}
CANARY_PANELS['T_NK'] = CANARY_PANELS['T_NK'] + ('CD45',)
CANARY_PANELS['B_Plasma'] = CANARY_PANELS['B_Plasma'] + ('CD45',)
CANARY_EXCLUSIVE_GROUP = {'Endocrine': 'Endocrine_islet', 'Exocrine': 'Exocrine_epi', 'Vascular': 'Endothelial',
                          'Stromal': 'Stromal', 'T_NK': 'Lymphoid', 'B_Plasma': 'Lymphoid', 'Myeloid': 'Myeloid',
                          'Neural': 'Neural'}
CANARY_PLACEMENT = {'T_NK': ('Imm-T-cell', ('CD3e',)), 'B_Plasma': ('Imm-B-cell', ('CD20',)),
                    'Myeloid': ('Imm-macrophage', ('CD68', 'CD163')),
                    'Exocrine': ('Exocrine', ('ECAD_uncompensated', 'panCK'))}
CANARY_MARGIN = 2          # Xenium DOMINANCE_MARGIN
CANARY_MIN_DEPTH = 2       # Xenium top_depth >= 2

# every label a cell no rule claims can carry: any canary panel marker or profile marker except ColIV
TOP_MARKER_LABELS = tuple(dict.fromkeys(
    top_marker_label(m) for m in list(TOP_MARKER_PROFILE) + [m for ms in CANARY_PANELS.values() for m in ms]
    if m not in TOP_MARKER_EXCLUDE))

# The hormone level, a departure (2026-10-17). Senior admits a cell inside an endocrine object when max(scaled INS,
# GCG, SST) clears the tissue p99 of that max. Two changes:
#  1. Each hormone is divided by its own ordinary-tissue p99 in the image (aligned values, so the image's UniFORM
#     factor cancels) before the max is taken, so the level no longer compares INS against GCG's background on the
#     cohort scale. Alone this changed little (Aab_6558 6,876 -> 6,888 endocrine cells).
#  2. Inside the objects the level is lowered, never raised, to where the cells clearing it are still at least 9:1
#     endocrine over background: the lowest level from which every bin up to the tissue p99 has a local false-
#     discovery rate <= LOCAL_FDR_MAX, background being the tissue distribution scaled to the object cells below
#     the tissue median (classify.object_hormone_level). The tissue p99 controls false positives among exocrine
#     cells at 1 %, which inside the objects is a global FDR of 0.1-0.2 % on every image -- far stricter than the
#     question needs, and on an image with a weak stain it cuts through the endocrine population: Aab_6558 (GADA+
#     only, C-peptide 8 ng/mL) has INS at 2.4x its tissue p99 in the objects against 15-86x on the other ND and
#     Aab+ images, a distinct INS-high / GCG-low cluster below the p99, and 35 % of its object cells Excluded.
#     With the local-FDR level well-stained images gain 0.2-3 % endocrine cells and the weak ones more (Aab_6558
#     22 %, T1D_6533 19 %, T1D_6534 9 %). The subtype winner is unchanged (highest scaled hormone).
HORMONE_LEVEL_DEPARTURE = ('max over hormones of aligned / own tissue p99; level = tissue p99 of that max, lowered '
                           'inside the objects to local FDR <= LOCAL_FDR_MAX')
LOCAL_FDR_MAX = 0.10
HORMONE_FOLD_BINS = np.exp(np.linspace(np.log(0.1), np.log(20.0), 41))   # log bins ~14 % wide

# Constants, unchanged from Senior.
REFERENCE_QUANTILE = 99.0
ENDOCRINE_QUANTILE = 1.0
MIN_REFERENCE_CELLS = 500
LARGE_ISLET_MIN_CELLS = 10
SIZE_STRATA = ('1', '2-3', '4-10', '>10')
IMMUNE_REQUIRE_MIN = 0.10      # the former single positivity level (v15 and before); refine.py's immune leaf still reads it
# Owner, 2026-10-28 (approved plan; output/positivity_method_comparison_20261028/): positivity is per image on the log
# scale, as the published pipelines do. A marker is positive when log1p(value) exceeds the image's background peak by
# POSITIVITY_K background widths (typing/positivity.py); the width is pooled over the images (the median) because
# REDSEA's clamp at 0 makes it unreadable in a third of them. The former level, scaled >= 0.10 on the linear cohort
# scale, sat at 3-170x the image background depending on the marker.
POSITIVITY_K = 3.0
POSITIVITY_WIDTH_POOLED = True
# Markers most cells carry have no background peak in tissue; their level is the empty-slide null (owner, 2026-10-28:
# "the background outside the tissue [helps] with the high expressers like ECAD").
MAJORITY_MARKERS = ('ECAD', 'ECAD_uncompensated', 'panCK', 'KRT8-18', 'BCatenin', 'BActin', 'epCAM', 'NaKATPase')
# How a majority marker's level is read from its own expressing population, floored at the empty-slide null
# (positivity.py): 'edge' = the population's lower edge, its peak minus 3 right-side widths (unimodal, E-cadherin);
# 'gate' = the two-component mixture gate on the nonzero log values (bimodal, panCK). Measured on v15's populations
# (output/positivity_method_comparison_20261028/epithelial_level_candidates*.csv): E-cadherin's edge keeps 97 % of
# the rule-typed Exocrine cells and 12 % of the rule-typed immune cells above it (the null: 99 % and 31 %), while
# panCK separates the two poorly at any level (an acinar cell's cytoplasm spills into its neighbor's cytoplasm
# compartment): null 87 % / 37 %, gate 52 % / 8 %.
MAJORITY_LEVEL = {'ECAD': 'edge', 'ECAD_uncompensated': 'edge', 'panCK': 'gate'}
# So E-cadherin alone is the epithelial negative of every immune rule (compensated, 'ECAD') and the Exocrine rule's
# marker (uncompensated, UNCOMPENSATED_COPIES); panCK keeps its level for the top-marker step and the positivity
# flags only.
NEGATIVES_NOT_ASKED = ('Pan_Cytokeratin',)
# Either measurement at its epithelial level makes a cell Exocrine: REDSEA errs one way on a junction (it zeroes a
# real one, never invents one), so a compensated positive stands and the uncompensated value adds the cells it zeroed.
# With the uncompensated value alone, 57,052 of v16's rule-typed Exocrine cells fell to the top-marker step.
EXOCRINE_RULE_MARKERS = ('ECAD', 'ECAD_uncompensated')
EMPTY_MARGIN_UM = 20.0
EMPTY_WINDOW_PX = 19           # 7.2 um square at 0.3774 um/px, about a cell's measurement area
EMPTY_WINDOWS = 20000
EMPTY_QUANTILE = 99.9
EMPTY_SEED = 0
SUPPORTING_LEVEL = 0.8715
MAJORITY_AREA = 0.5
CONFIRMED_MIN = 0.35            # islet_identity.REFERENCE_MIN
CONFIRMED_MARGIN = 2.0          # islet_identity.REFERENCE_MARGIN
UNDECIDED_MAX = 0.10
SHAPE_FLOOR = 1.0
ISLET_EXOCRINE_REASON = 'exocrine_inside_islet_border'
NO_MEASUREMENT_REASON = 'no_marker_values'
ENDOCRINE_OBJECTS = ('Islet', 'EndocrineCluster', 'EndocrineSingle')

# Senior's flattened immune requirements (celltype/rules.requirements over phenotype_workflow.csv, sha256
# b010a3f9...), dumped from Senior's own code on 2026-10-16: name -> (all_of, any_of groups, none_of).
SENIOR_IMMUNE_RULES = {
    'Imm-B-cell': (['CD20', 'CD79a'], [], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-NK-cell': (['CD56'], [['CD57', 'Granzyme_B'], ['CD11b', 'CD11c', 'CD20', 'CD209', 'CD66', 'CD68', 'CD79a', 'Iba1', 'MPO']], ['CD3e', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-T-cell': (['CD3e'], [], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-dendritic-cell': (['CD11c'], [['CD209', 'HLA_DR', 'IDO1']], ['CD68', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-macrophage': ([], [['CD68', 'Iba1'], ['CD11b', 'CD11c', 'CD20', 'CD209', 'CD66', 'CD68', 'CD79a', 'Iba1', 'MPO']], ['CD3e', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-neutrophil': (['CD66'], [['CD11b', 'MPO']], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Imm-plasma-cell': (['CD38', 'CD79a'], [], ['CD20', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Mac-M1-like': (['HLA_DR'], [['PD_L1', 'iNOS'], ['CD68', 'Iba1'], ['CD11b', 'CD11c', 'CD20', 'CD209', 'CD66', 'CD68', 'CD79a', 'Iba1', 'MPO']], ['CD163', 'CD3e', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'Mac-M2-like': (['CD163'], [['CD206', 'M2Gal3'], ['CD68', 'Iba1'], ['CD11b', 'CD11c', 'CD20', 'CD209', 'CD66', 'CD68', 'CD79a', 'Iba1', 'MPO']], ['CD3e', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-CD4-helper': (['CD3e', 'CD4'], [], ['CD8', 'E_cadherin', 'EpCAM', 'FOXP3', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-CD8': (['CD3e', 'CD8'], [], ['CD4', 'E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-CD8-cytotoxic': (['CD3e', 'CD8'], [['CD107a', 'Granzyme_B']], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-exhausted': (['CD3e', 'PD_1'], [['LAG3', 'TOX', 'VISTA'], ['CD4', 'CD8']], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-progenitor-memory': (['Bcl_2', 'CD3e', 'TCF_1'], [['CD4', 'CD8']], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
    'T-regulatory': (['CD3e', 'CD4', 'FOXP3'], [['CD39', 'ICOS']], ['E_cadherin', 'EpCAM', 'Ker8_18', 'Pan_Cytokeratin']),
}
# Rules this panel adds to Senior's (owner, 2026-10-16: "add the CD45 rule"). CD45 is on this panel and not on
# Senior's. It is a FALLBACK: asked only of cells no Senior definition claimed, before the supporting tier, with
# Senior's epithelial negatives -- so a CD3e+ or CD68+ cell keeps its Senior definition, and a leukocyte the panel
# cannot name further (Senior's B, NK, dendritic, neutrophil and plasma rules are unavailable here) is Immune
# rather than falling through to Exocrine.
ADDED_IMMUNE_RULES = {'Imm-CD45-leukocyte': (['CD45'], [], ['ECAD'])}
# The rule table's hierarchy for immune_subtype.py.
IMMUNE_LEVEL2 = ('Imm-B-cell', 'Imm-NK-cell', 'Imm-T-cell', 'Imm-dendritic-cell', 'Imm-macrophage',
                 'Imm-neutrophil', 'Imm-plasma-cell')
# Owner, 2026-10-27: no M1/M2-like labels; a macrophage is described by its marker profile instead (MACROPHAGE_PROFILE).
IMMUNE_LEVEL3 = {'Imm-T-cell': ('T-regulatory', 'T-exhausted', 'T-progenitor-memory', 'T-CD8-cytotoxic',
                                'T-CD8', 'T-CD4-helper')}
# Owner, 2026-10-27: "CD68+CD163+Macrophage_LGALS3+CD4+ with each plus also able to be a minus" (LGALS3 as Senior's
# galectin+ macrophage; CD4+ macrophages). Each marker is + at scaled >= IMMUNE_REQUIRE_MIN, else -. In text the
# plus and minus are superscripts; the data label keeps ASCII + and -.
MACROPHAGE_PROFILE = (('CD68', 'CD163'), ('LGALS3', 'CD4'))


def macrophage_label(positive: dict) -> str:
    """The macrophage label from {marker: bool}, e.g. 'CD68+CD163-Macrophage_LGALS3+CD4-'."""
    head, tail = MACROPHAGE_PROFILE
    sign = lambda m: m + ('+' if positive[m] else '-')
    return ''.join(map(sign, head)) + 'Macrophage_' + ''.join(map(sign, tail))


def macrophage_labels() -> list[str]:
    """Every label a macrophage can take (CD68 or CD163 is positive by definition), in a fixed order."""
    import itertools
    head, tail = MACROPHAGE_PROFILE
    out = []
    for h in itertools.product((True, False), repeat=len(head)):
        if not any(h):
            continue
        for t in itertools.product((True, False), repeat=len(tail)):
            out.append(macrophage_label(dict(zip(head + tail, h + t))))
    return out
# Owner, 2026-10-27: Senior's T-CD8 forbids CD4 and T-CD4-helper forbids CD8, so a T cell positive for both kept
# no subtype; and "once there is positivity for a broad defining marker, we do not need the 0.10 positivity level".
# A T cell (CD3e+) that neither rule names takes the subtype of the higher of the two markers, at any level, compared
# on the scaled values before the 0-1 clip (so cells clipped together at the floor or the top are still ordered).
T_SUBTYPE_BY_HIGHER = {'T-CD8': 'CD8a', 'T-CD4-helper': 'CD4'}
# Owner, 2026-10-27 ("I don't believe there is biological support for exact ties ... a new review of the total marker
# profile is indicated"): a T cell whose CD8a and CD4 are equal after the clip is reviewed on its whole profile. It stays
# a T cell when its highest marker among T_TIE_PROFILE_MARKERS (scaled, before the clip) is one of T_TIE_KEEP_MARKERS;
# otherwise the chain types it again without the T-cell rules (classify.py, resolved_by 'tie_profile_*').
T_TIE_KEEP_MARKERS = ('CD3e', 'CD8a', 'CD4', 'CD45')
# Refinement (typing/refine.py) defining markers, in this panel's names.
OWN_MARKERS = {
    'Nerve': ('B3TUBB',), 'Lymphatic': ('PDPN',), 'Blood vessel': ('CD31', 'CD34'),
    'Beta': ('INS',), 'Alpha': ('GCG',), 'Delta': ('SST',), HORMONE_NEGATIVE: ('CHGA',),
    'Immune': ('CD3e', 'CD68', 'CD20', 'CD45'), 'Exocrine': ('ECAD', 'panCK'),
}


def codex_rule(rule: tuple) -> tuple[tuple, dict]:
    """A Senior requirement in this panel's names, with every term this panel cannot ask removed.

    Returns (rule, removed). A rule whose required positives are all gone is returned with an empty
    all_of AND an empty any_of; `available_rule` marks such a rule unavailable rather than letting it fire on
    nothing. An OR group that loses every member is dropped, which also makes the rule unavailable when the
    group was its only positive evidence.
    """
    drop = set(NOT_ON_PANEL) | {s for s, c in SENIOR_TO_CODEX.items() if c in NOT_ON_EVERY_IMAGE}
    removed = {'all_of': [], 'any_of': [], 'none_of': []}
    all_of, any_of, none_of = rule
    keep_all = []
    for m in all_of:
        (removed['all_of'] if m in drop else keep_all).append(m)
    keep_any = []
    for group in any_of:
        kept = [m for m in group if m not in drop]
        removed['any_of'] += [m for m in group if m in drop]
        if kept:
            keep_any.append(kept)
        else:
            removed['any_of'].append('<whole group>')
    keep_none = []
    for m in none_of:
        (removed['none_of'] if m in drop or m in NEGATIVES_NOT_ASKED else keep_none).append(m)
    unavailable = bool(removed['all_of']) or '<whole group>' in removed['any_of']
    out = ([SENIOR_TO_CODEX[m] for m in keep_all], [[SENIOR_TO_CODEX[m] for m in g] for g in keep_any],
           [SENIOR_TO_CODEX[m] for m in keep_none])
    removed['unavailable'] = unavailable
    return out, removed


# Rules kept in a reduced form although a required positive is not on the panel, each with its reason. None.
ACCEPTED_REDUCTIONS: dict[str, str] = {}

# Rules asked in a substituted form on this panel, each with its reason (owner, 2026-10-17: "B cell should be" a
# cell type). Senior's Imm-B-cell needs CD20 and CD79a; CD79a is not on this panel. CD20 alone at scaled 0.10
# called 1,392-34,773 B cells per image -- CD20's background stays 11-fold apart between images after UniFORM --
# against 0-1,182 in the earlier typing. CD45 takes CD79a's place as the second required positive, and T-cell,
# macrophage and epithelial markers must be absent: 2-815 per image, the T1D donors highest (T1D_6533 785,
# T1D_6534 815, T1D_6520 394), close to the earlier typing's 0-1,182 (651, 1,182, 357). Reading CD20 over each
# image's own epithelial background (p99) instead left 0-10 on 11 of the original 15.
SUBSTITUTED_RULES: dict[str, tuple] = {
    'Imm-B-cell': ((['CD20', 'CD45'], [], ['CD3e', 'CD68', 'ECAD']),
                   'CD79a is not on this panel; CD45 is the second required positive, CD3e and CD68 must be '
                   'absent (owner, 2026-10-17)'),
    # Senior's Imm-macrophage needs CD68 or Iba1; Iba1 is not on this panel, which left CD68 alone. Owner, 2026-10-27:
    # a macrophage is CD68+ or CD163+ (so a CD68- CD163+ macrophage exists), with Senior's T-cell and epithelial
    # negatives. Senior's second group (any myeloid or B lineage marker) adds nothing beyond the first here.
    'Imm-macrophage': (([], [['CD68', 'CD163']], ['CD3e', 'ECAD']),
                       'Iba1 is not on this panel; a macrophage is CD68+ or CD163+, CD3e and the epithelial markers '
                       'absent (owner, 2026-10-27)'),
}

# Structure cell calls (owner, 2026-10-17: "Blood vessel, lymphatic, and nerve should not be cell types"). The
# Nerve, Capillary and Lymphatic objects stay -- the Microenvironment layer and each cell's in_<structure> /
# crossed fraction -- but no cell is typed from them: a cell inside one keeps its own identity through the
# endocrine, immune and Exocrine questions, as Senior's correction C already did for cells whose top marker was
# not the object's. False reverses that.
STRUCTURE_CELL_CALLS = False


def codex_rules() -> tuple[dict, dict]:
    """(available rules in this panel's names, per-rule record of what was removed and whether it survives).

    A rule that loses a REQUIRED positive (an all_of marker, or every member of an OR group) is unavailable:
    asking what remains would be a different question (Senior, immune_subtype.restrict_required: "reporting
    one as the other would let a donor that cannot be asked look like a donor that answered"). The exceptions
    are ACCEPTED_REDUCTIONS, which keep positive evidence on this panel.
    """
    rules, record = {}, {}
    for name, rule in SENIOR_IMMUNE_RULES.items():
        if name in SUBSTITUTED_RULES:
            rules[name], why = SUBSTITUTED_RULES[name]
            record[name] = {'substituted': True, 'rule_on_this_panel': rules[name], 'reason': why,
                            'senior_rule': rule, 'unavailable': False}
            continue
        r, removed = codex_rule(rule)
        if removed['unavailable'] and name in ACCEPTED_REDUCTIONS and (r[0] or r[1]):
            removed['unavailable'] = False
            removed['accepted_reduction'] = ACCEPTED_REDUCTIONS[name]
        record[name] = removed
        if not removed['unavailable']:
            rules[name] = r
    return rules, record


def cohort_rows() -> list[dict]:
    return list(csv.DictReader(COHORT.open()))


def manifest(out: Path) -> dict:
    return {e['image_name']: e for e in json.loads((out / 'project_manifest.json').read_text())['entries']}


def channels(out: Path, entry_id) -> list[str]:
    server = json.loads((out / 'qupath' / 'data' / str(entry_id) / 'server.json').read_text())
    return [c['name'] for c in server['metadata']['channels']]


def typing_dir(out: Path) -> Path:
    return out / 'typing'
