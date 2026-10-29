# Cell-typing rulings ledger

Every owner ruling on cell typing for the uniform 20-donor CODEX run, from this repository's sessions and the
Islet-Explorer-Senior sessions, verbatim with dates. **Read this before any typing work and before asking the owner
anything about typing**; ask only what no ruling here covers, and name the ruling you applied. A later ruling supersedes
an earlier one on the same point. Started 2026-10-28 with the final typing plan (`docs/final_typing_plan_20261028.txt`).
Sources: session transcripts of 2026-10-15 to 10-28 (this repository) and the Senior memory/transcripts; "[Senior]" marks
a ruling given in the Islet-Explorer-Senior project, which applies here as a principle.

Status: **in force**, **superseded** (by the ruling named), or **agent choice** (an agent's decision or an option the agent
wrote; not the owner's words, kept or changed on evidence).

## Measurement and positivity

| Date | Ruling (verbatim) | Status |
|---|---|---|
| 09-15 [Senior] | "We don't pool donors, the cell typing is per donor unless you can find evidence to the contrary." | in force (one image per donor here) |
| 10-19 [Senior] | "We will never make a comparison or do an analysis based on high or low levels of a marker; the vastly more important aspect is to have clear background determination to set binary positivity thresholds. Therefore actual max and min values don't matter between donors." | in force |
| 10-19 [Senior] | "the problem with your method is that one rule doesn't fit all situations." / "why would you take the densest point as the background peak? why are you treating all markers regardless of expression density the same?" | in force |
| 10-19 [Senior] | "Why are you using the typing's labels? This is supposed to be done before typing so the typing is more accurate" / "You need to use the values from after REDSEA." | in force |
| 10-26 [Senior] | "There can be no acceptance based on these outputs. Marker acceptance is tied to it's affects on object and cell type positivity." | in force |
| 10-16 | "For 5, there is a difference between actual background and variation between images." | in force |
| 10-18 | "We should be using the same processing tactics as the senior repo which utilizes REDSEA for spillover correction." | in force |
| 10-27 | "there should be no such thing as 'no cytoplasm'. This was already determined." / "You should use the whole-cell fallback" | superseded by 10-28 (one measurement per marker; whole cell exists for every cell) |
| 10-27 | Answer "Floor at zero" (an option the agent wrote) | superseded (v16) |
| 10-28 | "Why is it that almost every other pipeline uses normalization and scaling but you think it not necessary?" | in force (log scale, per image) |
| 10-28 | "THESE ARE CODEX IMAGES SO OF COURSE THEY ARE ALREADY BACKGROUND SUBTRACTED. ... I would expect the background outside the tissue to only help with the high expressers like ECAD" | in force |
| 10-28 | Asked to approve the v16 positivity method (mode + 3 pooled widths): "Discuss first" | v16/v17 method never approved |
| 10-28 | "We can't have two ways of typing. We need to use whole cell and/or membrane and/or nucleus for all measurements." | in force |
| 10-28 | "Don't use markers that are not in the panel for every image." | in force |
| 10-28 | "Also some images look like they need a despeckle filter prior to cell intensity measurements" | in force: every channel of every image through a 3x3 median (ImageJ Despeckle) before the cells and the empty-slide nulls are measured (`redsea_cells.py`); the filter choice, its measured effect and the evidence are the agent's (`output/typing_final_checks_20261028/despeckle/`) |
| 10-28 | "You have faster processing of information and access to more data so you can make the determinations." (about the positivity pass criteria) | in force: validation criteria are the agent's, recorded with sources and reported with results |
| 10-28 | "Don't flag the images or markers, flag the methods you are using so they will work with the images and markers" (on reading the Step 4 flags: Aab+/T1D failures listed, CD34/CK19/PGP9.5/VIM dropped, CD68 without a level on two images) | in force: a method that fails on an image or marker is revised, for every image, until it works on all of them; reports name the method and how it fails, not the image or marker. Supersedes the plan's "left out of typing in all 20 images" and the Step 4 "flagged for review" |

## Marker roles and method

| Date | Ruling (verbatim) | Status |
|---|---|---|
| 10-16 | "D3 and D4 and D5, use the cell typing from Islet-Explorer-Senior as much as possible" | in force, limited by 10-28 (Senior's lists were PhenoCycler-founded) |
| 10-17 [Senior] | "Refined is not a sensitivity analysis, it is part of the basic cell phenotyping: cells find their broad markers first to facilitate their true identity in a possible subclass" | in force (one hierarchical method) |
| 10-19 [Senior, relayed] | "There is no point in checking for impossible combinations in an unsupervised way when we have clear a priori knowledge." | in force |
| 10-28 | "The lineage and shared lists are not correct for the CODEX data because they were founded on the phenocycler data. The lineage and shared lists are not static and will be panel and tissue specific which needs to be clear for the pipeline as a whole. You have all the canonical marker definitions available to you in literature and general scientific knowledge so use it to correct the rules." | in force: `scripts/uniform/typing/registry/` |
| 10-28 | "The top marker needs to use the canary marker logic that the Xenium processing uses for canary genes" | in force, extended to every cell by the 10-28 request "to confirm cell identities" |
| 10-28 | "Lets do the updated processing and see what's left before making any other calls" | in force |
| 10-01 [Senior] | "use all non-diabetic/controls for testing classifier since the endocrine differences with disease interfere" | in force (shared settings set on the 7 ND images) |

## Classes and labels

| Date | Ruling (verbatim) | Status |
|---|---|---|
| 10-16 | "Cell quality check should drop cells without a nucleus, and that should take care of small cells, but do not drop irregular cells." | in force |
| 10-16 | "add the CD45 rule" | in force |
| 10-16 | "CHGA is supposed to be an aid to endocrine identification, not a blocker ... the GCG and SST can be used for typing" | in force |
| 10-16 | "A6, we need to drop cell classes that indicate neural, blood vessel or lymphatic and use the objects instead" | in force |
| 10-17 | "Blood vessel, lymphatic, and nerve should not be cell types, but B cell should be." | in force; the B-cell CD45 requirement is the agent's (canonical), not these words |
| 10-17 | "There should not be a phenotype listed that is a general class when there are subtypes" | in force (display) |
| 10-18 | "These hormones are mutually exclusive so there needs to be resolution or exclusion based on the approved islet cell typing methods. The 'only' and 'at all' designations are not biology, they are indicators of technical problems that need to be solved. We need to add a CHGA+INS-GCG-SST- cell type." | superseded 10-29 (highest scaled hormone, no exclusion); the 2x margin was an agent constant |
| 10-27 | Answer "Higher marker, no margin" to "How should a T cell positive for both CD8a and CD4 be named?" | in force for double positives |
| 10-27 | "Once there is positivity for a broad defining marker, we do not need the 0.10 positivity level." | premise (the 0.10 level) gone; T cells with neither CD8a nor CD4 positive: owner decision after the final run |
| 10-27 | "I don't believe there is biological support for exact ties, for these a new review of the total marker profile is indicated." / answer "Top marker, then unclipped" | moot (no ties with continuous values) |
| 10-27 | "We don't want to get into M1 and M2-like labels, but we can indicate macrophages and LGALS3 as defined in the senior workflow. We can also define CD4+ macrophages." / "Use CD68+CD163+Macrophage_LGALS3+CD4+ with each plus also able to be a minus" / answer "CD68 or CD163" | in force |
| 10-27 | "I don't want any more misclassified cells falling through the cracks. The established cell typing should classify cells appropriately and there should be no default of any class such as exocrine or immune. If typing has 'left over' cells and a check of top markers can't place them, call them by their top marker e.g. CD45+" | no cell-type default: in force; naming by top marker superseded by 10-28 "Get rid of naming by top marker logic" |
| 10-27 | "Any level and ColIV cannot be a cell type, if it is highest, it is a reason for exclusion." | "Any level" superseded (names only from positive markers, plan 10-28); ColIV part superseded by 10-28 |
| 10-27 | "Instead of Immune, the class should be CD45+" / "Keep them separate." | class name in force; "separate" moot in the hierarchy |
| 10-27 | Answer "Yes, leftovers only" (zero floor and top-marker step) | moot (the hierarchy applies to every cell) |
| 10-28 | "over 45 thousand ColIV cells is too many and likely another false result" | led to 10-28 below |
| 10-28 | "ColIV needs to be excluded from typing completely." | in force |
| 10-28 | "Each cell needs to start with a default class based on either tissue or islet parent." / "the object's class" | in force: start class `Tissue`, `Islet`, `EndocrineCluster` or `EndocrineSingle` |
| 10-28 | At the Step 4 stop: SMA by Barlow's known-negative 99th percentile, "1. use this"; B cells option (a) (CD20 above the 99.9th percentile of isolated epithelial cells, CD45+, single hematopoietic lineage, expected false calls reported), "2. use (a)" | in force (`docs/final_typing_handoff_20261028.md`) |
| 10-28 | "I thought the islets, clusters, and single endocrine cells were already accurately identified for all images." / "If the singles were already found and verified, don't undo that.  The endocrine typing didn't need to change" (during the final run's Step 6, on a draft that re-derived hormone names and moved single cells to other classes) | in force: the endocrine objects are used as given, and the established endocrine typing (v17's calls, unchanged method since v15) is kept decision for decision: Beta / Alpha / Delta, CHGA+ hormone-negative, unresolved subtype and its exclusions (hormone margin, exocrine inside the islet border, singles); the final run names no other cell endocrine (`identity.py`, labels file named_by `established`) |
| 10-28 | "I TOLD YOU TO USE WHOLE CELL FOR EVERYTHING AND I THOUGHT NO POOLED WAS A RULE" (on learning that the established endocrine typing compared cohort-scaled cytoplasm values) | in force: no cytoplasm and no pooled value anywhere; the endocrine typing is redone per image on whole-cell values inside the settled objects; whether the immune surface markers and E-cadherin may stay on the membrane (allowed by the 10-28 "whole cell and/or membrane and/or nucleus") is asked |
| 10-28 | "These are not compartment rules, they are lineage markers that we redid as part of this plan" (asked why no endothelial, stromal and other classes) | in force: CD31 names Endothelial, SMA Stromal, B3TUBB/GAP43 Neural cells (labels file); supersedes naming those cells by marker |
| 10-28 | "there should be no CD68+ or CD163+ cell classes since that would be covered under the macrophage typing" | in force with the next ruling: CD68 or CD163 positivity is the macrophage lineage, never a marker class |
| 10-28 | "YOU KNOW THERE IS NO SUCH THING AS A FUCKING CD8+CD163+ CELL!" (on a fallback that typed such a cell as a macrophage) | in force: a canary-level conflict between mutually exclusive compartments or lineages that the margin does not resolve gets no cell type; the cell is kept with the class "Conflict X/Y" ("Doublet suspected" for three or more groups); supersedes the plan's "keep the identity of the top compartment" |
| 10-29 | "The exocrine cell inside an islet border rule needs to be removed. It was designed to deal with inaccurate endocrine object borders which we fixed, so there should rarely if ever be a true exocrine marker profile inside object borders." | in force: an epithelial cell without an endocrine marker inside an islet is Exocrine (3,660 cells), never Excluded |
| 10-29 | "We need to be using the senior method of iteratively assigning cell types using the vessel, lymphatic, and nerve objects' checks for area and centroid. Blood vessel cells would be type "Endothelium Blood Vessel" and lymphatic "Endothelium Lymphatic" and "Neural"" | in force; supersedes 10-16 A6 and 10-17 ("Blood vessel, lymphatic, and nerve should not be cell types"): Senior's structure step (centroid in an object, more than half the area inside, the object's marker CD31 / B3TUBB / PDPN positive and the cell's top marker by fold over level, Blood vessel > Nerve > Lymphatic) before the canary conflicts; clusters and singles keep their object |
| 10-29 | "The 2x runner up still isn't working, see Object ID ddb348e9-c4bb-4ce0-9dc0-7a2e22f44533 in ND_6232" (GCG 2,054, SST 118, excluded under the margin) | in force: the 2x margin is read on Senior's typing scale per image, (x - level) / (p99.9 - level), floored at 0, unclipped above (10-27 "then unclipped"), not on folds over the levels; ddb348e9 Alpha, af2ee599 Delta; margin exclusions 11,656 -> 7,533 |
| 10-29 | "name islet cells by highest scaled hormone as senior's method does" | in force; supersedes the 10-18 resolution-or-exclusion and the 2x margin: the unique highest hormone among those above their levels, on the per-image typing scale, names the cell (Senior subpanel step 5); an exact tie keeps the start class with status tied_highest_hormones (none on the 20 images); the 7,531 margin exclusions become Alpha 2,977, Beta 2,525, Delta 2,029; HORMONE_MARGIN removed from the registry |
| 10-29 | "The canary and related checks are not working as expected. For Aab_6450 Object ID 1b2c9da1-000c-4730-9a3e-88e5ec21994c is listed as conflicting, but has very strong CD3e and CD45 which should have pushed it to T Cell CD8a"; "There are two markers for a T Cell CD8a and both are brighter and cover more cell area than any myeloid marker so this should not be a conflicting call." | in force (method repair): between the lineages of one compartment the canary margin is CANARY_LINEAGE_MARGIN 1, not 2 (compartments keep 2). Pan, shared and multi-lineage markers add to every lineage, so only 1-2 lineage-specific markers separate T, B and Myeloid; at 2 any canary-positive marker of another lineage made a conflict. 1b2c9da1 (T: CD3e, CD8a; Myeloid: CD163) is CD8a Tcell; 4,288 conflicts named on 20 images (CD8a Tcell 2,592, macrophage profiles 1,223); equal depth stays a conflict |
| 10-28 | "Get rid of naming by top marker logic" | in force: no cell is named by a marker ('<marker>+'); a cell its compartment and lineage do not name keeps its start class with a status (no_lineage_named, <lineage>_marker_only, shared_marker_only, no_marker_above_background); supersedes the 10-27 "call them by their top marker e.g. CD45+" and the plan's top-marker step. The pan-leukocyte label CD45+ (10-27 "Instead of Immune, the class should be CD45+") is a lineage label and stays |
| 09-18 [Senior] | Exocrine is one label (the Senior panel had no acinar or ductal marker) | in force for this run; premise differs here (CK19 on this panel): split is an owner decision after the final run |

## Display

| Date | Ruling | Status |
|---|---|---|
| 10-20 | "Don't use the word 'share' where it means proportion." | in force |
| 10-27 | Broad name first, superscript ±, "T Cell, B Cell", "Endocrine first"; macrophage marker checkboxes | in force |
| 10-26 | American English spelling | in force |

## Agent choices recorded as such (not the owner's words)

pooled background width (dropped 10-28); E-cadherin read twice (dropped 10-28); panCK dropped as an immune negative (v16;
superseded by the canonical registry); local-FDR hormone level (10-17; one hormone definition chosen by the endocrine-object
check); hormone 2x margin (withdrawn 10-29); Senior's supporting tier (removed 10-27, never approved); CD20 with CD45 for B cells (canonical);
K = 3 background widths; the canary detection level; the positivity pass criteria (delegated 10-28); the 3x3 median despeckle (10-28, on the owner's request); CD8a and CD4 known-negative levels (10-28 method repair after the despeckle: every background-peak method called touching CD8a+ structures in 5.1-6.1 % of ND_6479's isolated hormone cells; CD4 moved with CD8a so the higher-marker comparison reads both over one kind of level).
