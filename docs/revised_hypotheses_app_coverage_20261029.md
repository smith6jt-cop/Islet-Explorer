# Revised hypotheses using the full Islet Explorer dataset

29 October 2026. Current conceptual revision of the [27 October proposal](revised_hypotheses_beyond_barlow_20261027.md). R1–R6 retain their subjects; R7–R11 add exocrine tissue, lymphoid organization, Ki67 positivity, integration across app representations, and clinical relevance. This revision specifies questions and measurements; it does not report results or start analysis on the data being repaired elsewhere.

## Central scientific question

**Among endocrine objects with comparable composition and surrounding inflammation, which cellular identities and tissue configurations distinguish immune localization outside the object from localization within it?** How do those relationships vary across ND, Aab-positive and recent-onset T1D donors, across object sizes, and across locations within a section?

The study uses cells, endocrine objects, extracellular structures, physical coordinates and donor metadata together. Each has a defined role. A measurement can define identity, test a hypothesis, provide a denominator, or check an alternative explanation; every marker does not need its own mechanistic claim. Likewise, several plots of the same measurement are complementary views, not independent confirmation.

Barlow already examined endocrine and immune composition, core/peri trajectories, vascular and neural abundance, regional organization and lymphoid neighborhoods. The proposed contribution is the **relationship between specific macrophage profiles, measured interfaces and the surrounding epithelial, neural and vascular setting**, especially among otherwise comparable endocrine objects and among small objects outside the embedding. Barlow's acquisition of a marker does not mean every question involving it has been exhausted, and a different analysis alone does not establish novelty. [Barlow et al.](https://pubmed.ncbi.nlm.nih.gov/40232951/)

## Data and measurement commitments

- The intended cohort remains 20 donors: 7 ND, 7 Aab-positive and 6 T1D. Donor 78 is excluded. Donors, not cells or islets, are the biological replicates. Disease groups are compared directly; Aab positivity and pseudotime do not establish a temporal sequence or a clinical stage.
- Use the corrected, accepted root when it becomes available. On this date, `uniform-20donor-20261029-v18/uniform_release.json` says `assembled_pending_owner_review`, and further tissue and typing work is underway. Its schema informs this plan; its cell/object counts are not adopted as validated results. Historical v10 counts, coordinates and endpoint tables cannot be joined to changed identities by assumption.
- Follow the latest [typing rulings](typing_rulings_ledger.md): accepted exclusive identities, per-image marker positivity, retained unresolved/conflict records and structure objects. Do not reconstruct identity from an old generic gate, call CD163 a universal macrophage-specific marker, or assign M1/M2 labels.
- Use counts, presence, proportions, area, distances, geometry and accepted binary marker states. Cross-donor high/low marker-intensity comparisons are outside this plan. The app's scaled marker-expression display is available for inspection but is not a calibrated protein-abundance endpoint. This applies to INS and to the broader standing marker-measurement instruction.
- Treat blood vessels, lymphatics, nerves and Collagen IV as structure objects. Endothelial, Stromal and Neural **cell identities**, where present in the accepted typing, remain separate measurements. Membership or overlap with a structure is not itself a cell identity. ColIV is neither a cell identity nor a cellular positivity endpoint.
- Preserve the full tissue and endocrine-object census. UMAP/DPT eligibility is at least ten total core cells and does not define the whole biological population. Islets and Clusters remain distinct annotation classes. Singles contribute counts and counts per tissue area only, under the existing scope.
- Use the accepted core and tissue-clipped 20-µm peri band, plus explicitly defined outer tissue. Overlapping bands can contribute to separate object summaries but must be deduplicated for tissue-wide or donor cell totals. Missing, unmeasured and ineligible observations are not biological zeros.

## The revised hypotheses

Each hypothesis has a direct ND/Aab-positive/T1D comparison and a conditional comparison that addresses its biological question. Display both. A group difference does not establish a within-group relationship, and different within-group p-values do not establish a group difference.

### R1 Macrophage phenotype and position distinguish immune access

**Hypothesis.** At comparable surrounding T-cell and total macrophage abundance, LGALS3-positive macrophage profiles are more concentrated at interfaces with core-localized T cells than at interfaces where T cells remain outside. CD68/CD163 background, CD4 and HLA-DR positivity further distinguish that association.

**Test.** The main response is the proportion of local T cells in the core among T cells in core plus peri, restricted to objects with at least one local T cell. Show core and peri T-cell counts and area densities alongside this allocation. The primary macrophage contrast is LGALS3-positive versus negative within macrophages; CD4, the three CD68/CD163 backgrounds and HLA-DR are prespecified refinements. Separate CD3-defined CD4 and CD8 T cells where supported. Zero-exposure objects remain in the census and presence analyses rather than acquiring an invented allocation of zero.

**App use.** Phenotype composition/enrichment, the immune gradient using each donor's own tissue baseline, signed-distance profiles, Plot against another measure and donor distributions. UMAP/Leiden locate the same objects in endocrine composition space. DPT/PAGA show their position on the existing core/peri composition graph; correlations with a macrophage or T-cell feature used to build that graph are descriptive.

**Discriminating outcome.** A profile association beyond total burden would show information lost by pooling macrophages. Adequately precise absence would favor burden or measured anatomy for this endpoint. Neither result establishes antigen presentation, killing, clearance or protection. Direct pancreatic evidence supports a macrophage interpretation of many CD163-positive cells; DC3 presence in our nonmalignant samples remains unestablished. [Non-diseased pancreas study](https://pubmed.ncbi.nlm.nih.gov/31851923/); [pancreatic cancer DC3 study](https://pubmed.ncbi.nlm.nih.gov/37938368/).

### R2 Collagen IV configuration distinguishes local access from general inflammation

**Hypothesis.** Immune localization within endocrine objects is associated with the placement of ColIV discontinuities relative to vascular intersections and macrophage profiles, beyond total ColIV area and the amount of inflammation outside the object.

**Test.** Start with accepted ColIV object area in core, peri and outer tissue. Derive boundary coverage and discontinuities separately: area occupancy is not boundary coverage. Compare immune localization beside covered boundary segments, gaps near blood-vessel intersections and gaps away from those intersections, with segment length, surrounding immune cells and object size accounted for. Include the whole-tissue ColIV mesh so an endocrine-interface association can be distinguished from a general difference in matrix abundance.

**App use.** Structure-area versus diameter; disease and donor distributions; tissue maps; marker-independent phenotype distance profiles; selected-object segmentation and Viewer overlays. Plot the derived interface quantities on the fixed composition map and along DPT as external annotations. Segment-wise gap measurements require a derived geometry table and are not currently an app card.

**Discriminating outcome.** A site-specific relationship distinguishes interface configuration from overall matrix content. Equivalent behavior near all gaps weakens vascular specificity. A staining gap in a 2D section is not sufficient evidence of physical barrier failure. Basement-membrane loss at infiltration sites is already documented; that broad claim is not the new hypothesis. [Korpos et al.](https://pubmed.ncbi.nlm.nih.gov/23139348/)

### R3 Neural configuration distinguishes immune niches beyond nerve abundance

**Hypothesis.** GAP43-positive portions of accepted neural structures have different macrophage and T-cell proximity than GAP43-negative portions with comparable neural extent and vascular context.

**Test.** Use B3TUBB and accepted nerve geometry to define the neural setting; compare immune counts in specified surrounding regions and distances to the two neural configurations. Adjust or match for neural extent, local immune abundance and blood-vessel proximity. Report nerve area per tissue area alongside nerve area per endocrine area. PGP9.5 provides supporting localization but cannot independently define nerves because endocrine cells can express it.

**App use.** Neural structure area in core/peri, Plot against diameter and other measures, GAP43-positive distance profiles, tissue maps and image inspection. UMAP and the multi-feature DPT heatmap show whether these configurations vary among endocrine composition states. Nerve-segment-specific immune distances require derivation from existing geometry.

**Discriminating outcome.** Persistence after accounting for nerve abundance supports a configuration-specific association. Loss of the association after that adjustment favors abundance or geometry. GAP43 positivity is not evidence of new innervation, regeneration or functional neural activity. Nerve and capillary persistence despite beta-cell loss has already been studied. [Richardson et al.](https://pubmed.ncbi.nlm.nih.gov/36696598/)

### R4 Blood vessels lymphatics and surrounding stroma define different immune settings

**Hypothesis.** Blood-vessel-associated and lymphatic-associated neighborhoods differ in the outside/inside distribution of T cells and macrophage profiles, and SMA-positive stromal context modifies those relationships.

**Test.** Compare accepted blood-vessel and lymphatic object area and proximity to endocrine objects at comparable immune exposure and endocrine size/composition. Separate overlap with either structure from proximity to both. Within these settings, distinguish accepted Stromal/SMA-positive cells adjacent to vessels from those elsewhere, and retain CD34, VIM and PDPN positivity as contextual measurements within their accepted lineages. Do not classify every PDPN-positive cell as lymphatic endothelium or every SMA-positive cell as an arteriole or pericyte.

**App use.** Core/peri/core-plus-peri structure-area plots, cell-population counts and densities, donor curves, phenotype enrichment, tissue scatter and Viewer. Display blood-vessel, lymphatic, stromal and immune measurements together on the existing DPT heatmap; geometry not used in the fit is an external annotation.

**Discriminating outcome.** Different configurations at comparable structure amount would show that blood and lymphatic context contain distinct information. Similar configurations favor a common inflammatory setting. Static sections cannot establish ingress, egress, drainage or perfusion. Human peri-islet lymphatic localization is already described. [Human lymphatic study](https://pubmed.ncbi.nlm.nih.gov/26822093/)

### R5 Endocrine arrangement changes beta cell exposure beyond composition

**Hypothesis.** The observed arrangement of alpha, beta and delta cells changes the frequency of beta–immune near contacts relative to what is expected from the same endocrine composition and geometry.

**Test.** Use the already approved fixed-geometry, fixed-composition endocrine-label permutations. Hold immune positions and the spatial contact/near-contact definition fixed. Compare observed beta–immune adjacency with its permutation expectation; separately examine alpha–beta and delta–beta adjacency and beta-cell location at the object boundary. Preserve the CHGA-positive hormone-negative class rather than assigning those cells a former beta identity. Use contact opportunities and available immune-cell counts as denominators or adjustments.

This sharpens the earlier R5: simply showing that more perimeter touching alpha cells leaves less perimeter for something else would be partly a mathematical consequence. The revised endpoint tests whether the actual endocrine arrangement changes immune contact relative to the explicit spatial reference.

**App use.** Composition versus diameter, paired feature plots, saved composition UMAP/Leiden, linked single-cell maps and composition bars, and Viewer checks of selected interfaces. Adjacency and permutation residuals are derived research measurements; the app does not calculate them from a button.

**Discriminating outcome.** Observed contact differing from the fixed-composition reference supports information in arrangement beyond abundance. Similarity to the reference weakens that claim. A two-dimensional near contact does not demonstrate a physical synapse, endocrine paracrine signaling or prior protection. Broad architectural differences are already described. [Van der Heide et al.](https://pubmed.ncbi.nlm.nih.gov/41673026/)

### R6 Small endocrine clusters have distinct relationships to their local setting

**Hypothesis.** Small beta-bearing Clusters and larger Islets differ in how their tissue-area densities relate to macrophage profiles, matrix configuration and neural/vascular context.

**Test.** Use the full geometric census, including objects outside the embedding. Count objects in tissue regions with and without endocrine objects, using measured tissue area. Compare the size-by-setting relationship directly and retain separate beta presence, endocrine composition and object-size measurements. Do not equate a 2–9-endocrine-cell group with Barlow's fewer-than-ten-total-cell exclusion. Singles remain a separate donor count/count-density series; they receive no object-level neighborhood or trajectory model.

**App use.** Endocrine objects per donor; counts and counts/mm²; against-diameter plots and density contours; cumulative distributions; one curve per donor; islet/cluster filters; tissue maps; both nearest-Islet and nearest-Islet-or-Cluster distance views. The latter two answer different questions because small dispersed Clusters change which boundary is nearest.

**Discriminating outcome.** A size-specific association supports local context as an additional descriptor of small-object distribution. A common relationship across sizes favors a broader tissue effect. Scarcity alone cannot distinguish prior destruction, development, sectioning or detection. Small beta-rich object scarcity has already been reported. [Murrall et al.](https://pubmed.ncbi.nlm.nih.gov/41223263/)

### R7 Epithelial source distinguishes extraislet immune associations

**Hypothesis.** LGALS3- and HLA-DR-positive epithelial regions, particularly CK19-positive duct-associated regions, have different nearby T-cell and macrophage configurations from other exocrine regions at comparable local immune abundance.

**Test.** Within accepted Exocrine/epithelial identities, compare CK19-positive and negative cells and their local setting, then LGALS3 and HLA-DR positivity within those strata. The primary comparison is local CD4 T-cell abundance/proximity around HLA-DR-positive versus negative CK19-positive epithelial cells, accounting for ductal extent and surrounding immune burden. Compare the epithelial LGALS3 association with the macrophage-source association from R1. Relate nearby endocrine immune allocation as a secondary endpoint. ECAD, panCK and BCatenin support epithelial localization; CD44, KRT14, NaKATPase and BActin provide secondary binary-state and spatial-coverage profiles rather than identities or presumed functional states.

**App use.** Tissue-wide marker-positive distance/composition profiles and Viewer channels establish where the measured signal lies. Plot and donor views show exported lineage-specific positive fractions. Cell identity × marker positivity × local proximity needs a derived join; the existing marker-positive profile pools cell types and cannot by itself establish the cellular source. CK19 is a measured phenotype refinement here, not an unauthorized change to the Exocrine typing label.

**Discriminating outcome.** Different associations for epithelial and macrophage LGALS3/HLA-DR would show why pooled marker positivity is insufficient. Similar associations after accounting for nearby immune cells weaken source specificity. Epithelial HLA class II and even cytokine-induced antigen presentation have already been studied in human ductal cells; this CODEX panel tests spatial relationships, not antigen presentation. [Human ductal-cell study](https://pubmed.ncbi.nlm.nih.gov/42096071/)

### R8 B cell organization adds information beyond B cell abundance

**Hypothesis.** At similar local B- and T-cell numbers, B-cell-associated neighborhoods differ in T-cell proximity, macrophage profile and CD35-positive organization, and these differences associate with beta-containing endocrine objects.

**Test.** Separate B-cell presence and density from B-to-T spatial proximity and concentration around endocrine objects. Report CD35 positivity separately within B cells, macrophages and unresolved/nonimmune regions; use the images to distinguish a reticular staining pattern from positivity on individual leukocytes. Compare B-associated configurations with T/macrophage-dominant configurations and report which donors supply observations. Age and autoantibody number are secondary donor-level comparisons, not grounds for discovering many small endotypes.

**App use.** B-cell core/peri proportions and own-tissue enrichment; donor heatmaps; the immune gradient and signed-distance profiles; tissue scatter; the immune-versus-beta-content card; CD35 image overlays. B-to-T proximity and reticular morphology require explicit derived measurements or image annotation. A selector that suppresses very rare measures must not be interpreted as evidence that the population is absent.

**Discriminating outcome.** An association after accounting for B-cell abundance supports organization as additional information. Loss after abundance adjustment favors a cell-number explanation. Call these B-cell-associated aggregates or neighborhoods: CD35 alone establishes neither follicular dendritic-cell identity nor a mature tertiary lymphoid structure. Human pancreatic tertiary lymphoid organs are already documented, as are Barlow's lymphoid neighborhoods. [Human TLO study](https://pubmed.ncbi.nlm.nih.gov/33912981/)

### R9 Ki67 positivity is linked to lineage and local context

**Hypothesis.** The proportion of Ki67-positive cells within a lineage differs between endocrine-associated and outer-tissue settings, with different patterns across disease groups.

**Test.** The first contrast is macrophage Ki67-positive fraction in the peri region versus outer tissue, summarized within donors; T cells, B cells, endocrine subtypes, Exocrine, Endothelial and Stromal cells are secondary lineage-specific contrasts. Each denominator is the number of callable cells of that lineage in that region. Report positive counts and donors as well as fractions. Compare endocrine Ki67-positive events separately in beta, alpha, delta and CHGA-positive hormone-negative cells, with age shown explicitly. Sparse events receive descriptive estimates rather than an unstable multi-covariate model.

**App use.** Linked cell maps and Viewer support localization; existing donor, disease, diameter and clinical-age plot forms organize exported Ki67 summaries. Ki67 is in the common upstream panel/registry, but the inspected v18 `cell_distance.csv` has no Ki67 positivity column. This part needs the accepted upstream Ki67 calls joined to the analysis tables; it is not already available as a complete app distance analysis. The inspected manifest records whole-cell Ki67 measurement, so this plan does not describe it as a validated nuclear-only call.

**Discriminating outcome.** A regional difference within the same lineage supports niche-associated Ki67 positivity. A pooled difference that disappears within lineages reflects changing composition. Ki67 positivity is not a measured proliferation rate or proof of beta-cell replacement. Human studies have already examined islet Ki67 positivity, including hormone-negative proliferating cells; our panel cannot assign those cells the ARX/SOX9 identities reported elsewhere. [Recent-onset study](https://pubmed.ncbi.nlm.nih.gov/20532863/); [alpha-related hormone-negative cells](https://pubmed.ncbi.nlm.nih.gov/29326366/).

### R10 Tissue context distinguishes objects within the same composition state

**Hypothesis.** Disease-associated differences involve the mixture of immune and structural configurations within endocrine composition states, as well as the frequency of the states themselves.

**Test.** Use the fixed accepted composition UMAP and a prespecified Leiden resolution to describe states, then compare R1–R9 measurements within states or among objects with comparable defining composition. Show state composition, the number and identity of contributing donors, donor-specific state occupancy, and physical locations in tissue. Inspect neighboring objects with similar composition but different immune allocation and quantify those comparisons with the established spatial measurements. Unannotated tissue regions are not automatically anatomical lobules.

Use the corresponding app DPT and PAGA as a second descriptive representation of core/peri composition. Overlay matrix, vessel and neural geometry and other measurements that are absent from the fitted feature matrix. Record the actual feature list because corrected labels can change the graph inputs. A macrophage profile, immune fraction or cycling label that enters the fit is not an independent validation of that trajectory. The beta-ratio orientation is a modeling convention. Do not fit a new macrophage-defined trajectory as part of this revision.

**App use.** All UMAP colorings; Leiden composition, donor-status and donor membership views; raw and rank DPT; overall and disease-specific LOESS; the 25-bin status/feature strip; multi-feature heatmaps; PAGA; tissue maps; linked drill-down. Use common coordinates, scales, bins and fitted cohort across comparison panels. Rank DPT is not recalculated after choosing a marker or a disease facet.

**Discriminating outcome.** Differences within a composition state support information in tissue context beyond changing state frequency. Differences confined to state occupancy favor composition as the main organizing description. This is the integration of R1–R9, not another independent discovery test on the same endpoints. Low-dimensional visual separation, graph branches and LOESS bands do not establish disease trajectories or donor-level statistical significance.

### R11 Spatial context contributes to clinical associations beyond sampled beta quantity

**Hypothesis.** A prespecified donor summary of local immune organization is associated with recorded C-peptide beyond sampled beta-cell quantity and the selected clinical context.

**Test.** The first candidate is the donor summary of R1's LGALS3-positive macrophage enrichment near endocrine objects relative to that donor's own outer tissue. Compare a baseline model containing sampled beta quantity and the previously approved limited clinical covariates with the same model plus this one spatial summary. Reuse the approved nested log-C-peptide, donor-held-out prediction-error and censored-regression approaches; the new spatial predictor is a declared change from the historical fixed score. Fix its definition before seeing its clinical association. Recalculate every fitted transformation or feature-selection step within training donors if any such step is used.

Show C-peptide censoring explicitly; an app point plotted at a detection limit is not an uncensored assay measurement. HbA1c and admission/peak glucose describe glycemic context. Age, sex and body size are limited sensitivity variables. Autoantibody count/specificity characterize Aab-positive donors; T1D duration has little variation in this recent-onset cohort. Donation/recovery variables assess technical or terminal-illness alternatives. Do not fit all available donor fields simultaneously to 20 observations.

**App use.** Against a donor clinical value, donor medians and quartiles, clinical/demographic UMAP coloring, cohort metadata and donor-level summaries. Objects displayed behind each donor's clinical value remain repeated tissue observations, not additional clinical subjects. Clinical modeling and censoring correction are exported research analyses, not provided by the Plot card's Spearman correlation.

**Discriminating outcome.** An association and improved held-out error would support incremental descriptive value in this cohort, with wide uncertainty expected. Lack of improvement would favor sampled beta quantity and selected clinical context. Neither establishes prognosis, future progression or a clinically validated biomarker. This is exploratory clinical integration.

## Coverage of every common marker

The registry lists 32 markers present on every image. DAPI is an additional image/segmentation reference. The inspected distance export has 30 marker-positivity columns: Ki67 is absent and ColIV is handled as structure objects. That export and the upstream panel are not interchangeable. The following mapping gives each common marker a scientific or measurement role without inventing unavailable functions.

| Marker or markers | Main use | Hypotheses and limits |
|---|---|---|
| INS, GCG, SST | Exclusive beta, alpha and delta identity; counts, presence, proportions and arrangement | R1–R6, R9–R11. No quantitative hormone-content endpoint. |
| CHGA | Endocrine support and the CHGA-positive hormone-negative population | R5, R6, R9, R10. Not evidence of former beta identity; PPY and lineage-tracing markers are unavailable. |
| CD45 | Named immune populations and separately retained Hematopoietic CD45-positive/negative residual categories | R1, R7–R10. Keep the accepted immune aggregation definition explicit; residual/conflict categories are not silently added to T, B or macrophage totals. |
| CD3e, CD8a | Accepted T-cell identities and T-cell positioning | R1–R10. CD8a alone is not a new T-cell identity or proof of cytotoxic activity. |
| CD4 | CD3-defined CD4 T cells and CD4-positive macrophage profiles separately | R1, R7–R9. A shared marker does not establish a shared lineage. |
| CD20 | Accepted B-cell identities, density and organization | R8, with R1/R9/R10 context. Retain the current CD45/lineage requirements. |
| CD68, CD163 | Exclusive macrophage profiles | R1–R4, R7–R11. No M1/M2, resident/recruited or DC3 assignment from these markers. |
| LGALS3 | Binary macrophage-profile contrast and epithelial-source comparison | R1, R7, R11. Spatial positivity does not establish secretion or damaging/protective function. |
| HLADR | HLA-DR-positive fractions within accepted myeloid, B/T, epithelial and endothelial populations | R1, R4, R7, R8. Not a DC gate or an antigen-presentation assay. |
| CD35 | Cell-associated positivity and image-visible organization in B-cell-associated regions | R8. Reticular staining alone does not prove FDC identity or a mature TLO. |
| ColIV | Accepted tissue-wide and endocrine-associated matrix objects; area and derived boundary geometry | R2, R4, R6, R7, R10. No cell-positive fraction; old cellular ColIV flags are superseded. |
| CD31, CD34 | Blood-vessel/endothelial context; CD34-positive fractions within defined lineages | R2, R4, R6, R10. CD34 positivity can also be stromal; objects and cells have separate denominators. |
| PDPN | Accepted lymphatic geometry and positivity in defined endothelial/stromal cells | R4, R8, R10. PDPN-positive does not automatically mean lymphatic. |
| SMA, VIM | Stromal context, vessel-associated versus other locations, lineage-specific positivity | R2, R4, R6, R7, R10. Neither establishes fibrosis, a myofibroblast program or a pericyte identity alone. |
| B3TUBB, GAP43 | Accepted nerve geometry and GAP43-associated neural configurations | R3, R6, R10. No inferred nerve growth or function. |
| PGP9.5 | Supporting neural/endocrine localization | R3, R5. Interpret with structure and endocrine identity. |
| ECAD, panCK, CK19, BCatenin | Epithelial localization and CK19-positive refinement within accepted Exocrine cells | R5, R7. A whole-cell positive call does not measure junction integrity or Wnt activity. |
| CD44, KRT14 | Secondary epithelial/stromal state and spatial-coverage comparisons using accepted display-only positivity | R7, R10. No new identity, stemness, regeneration or EMT claim. |
| NaKATPase, BActin | Secondary lineage-specific positivity/coverage profiles and checks that apparent state differences are not broad shifts affecting many channels | R7, R10 and sensitivity analyses. Do not assume constant expression or use them as unvalidated normalizers; positivity is not pump activity or cytoskeletal function. |
| Ki67 | Positivity within named lineages, including endocrine hormone-negative cells | R9. Use the accepted upstream measurement; source-level availability does not mean the app already exports it. |
| DAPI | Nucleus support, cell and object geometry validation, image localization | All. A measurement-control role, not an independent biological hypothesis. |

KRT8-18, epCAM and Ghrelin are not common to every image and are excluded. No hypothesis depends on absent DC-lineage, exhaustion, Treg, IDO, cytokine or antigen-specificity markers. Historical protein embeddings remain provenance records rather than a substitute for the current composition representation. This CODEX plan does not assume matched Xenium, RNA, genotype or longitudinal measurements that the app has not supplied.

## Coverage of other data layers

| Data layer | Required role |
|---|---|
| Full-section cell identities and binary marker states | R1, R7–R9; establish cell source and the outer-tissue reference. Retain callable-cell counts and unresolved/conflict records. |
| Whole-cell and nuclear geometry, centroids and stable IDs | Image support, boundary/contact measurements and correct joins. Centroid proximity is not membrane contact. |
| Islet, Cluster and Single census; area and diameter | R5/R6; preserve nonembedded objects and zero-object tissue areas. Diameter is the area-equivalent diameter, not an independently measured anatomical diameter. |
| Core, 20-µm peri and outer tissue | Separate local access, exposure and the broader tissue baseline across R1–R9. |
| Structure objects and union area | R2–R4/R6/R10; intersect unions with the stated region, counting overlapping area once. Fragment count is not vessel/nerve number. |
| Structure-qualified cells, overlap/containment and passthrough fields | Where retained in the accepted root, describe which named cells occupy structures; keep these distinct from structure area and identity. Availability is recorded, not inferred from an old pipeline. |
| Signed distance to nearest Islet and nearest Islet or Cluster | R1/R6–R8; compare definitions explicitly. Saved object-center `min_dist_*` fields are not boundary distances. |
| Composition UMAP and Leiden labels | R5/R10 and common annotation of the other hypotheses; retain fitting features, resolution, eligible objects and donor contributions. |
| DPT, rank DPT, PAGA graph/layout and feature tables | R10 and external annotations for R1–R9; retain input features, root, orientation and software versions. No temporal or independent-validation claim from fitted inputs. |
| Raw multiplexed images and segmentation overlays | Assess the actual cells/structures underlying each selected example. Use all-image acceptance evidence from the repair work; do not replace it with attractive examples selected from significant results. |
| Clinical and demographic records | R11 plus limited prespecified confounding/sensitivity checks below. |
| Edge eligibility, tissue outlines, lineage confidence and missingness | Define measurement support and test whether conclusions depend on ambiguous or incompletely observed regions. This is substantive use of these data, not a new cell class or biological explanation. |

All donor fields are accounted for: disease group, autoantibody number and GADA/IA2A/mIAA/ZnT8A flags characterize disease context; C-peptide with censoring, HbA1c and admission/peak glucose characterize recorded metabolic context; age/age group, sex, height, weight, BMI/BMI percentile support demographic description and limited adjustment; pancreas weight describes organ context without converting sampled sections into whole-organ beta mass. ICU days, transit time, RIN, recovery type and cause of death support technical/terminal-illness sensitivity checks. Race/ethnicity and ABO group remain descriptive metadata unless a separate scientifically motivated comparison is specified; race/ethnicity is not a genetic mechanism. Donor IDs provide linkage and donor-contribution checks, never biological predictors. Preserve missingness; HbA1c recorded as zero is not a physiological zero, and duration is applicable only to T1D.

## Coverage of every app analysis and display

The current source defines four default tabs: Plot, Trajectory, Spatial and Viewer. Statistics and Inflammation are implemented but disabled by default. The tables below distinguish current app views from derived research measurements; they do not turn those modules on or change the app.

### Plot tab

Use both data layers: **Cell Populations** for counts, proportions and counts per area, and **Microenvironment** for structure area as a proportion of the stated region. Pair proportions with counts/densities because a falling denominator can raise a proportion without an increase in the numerator.

| Plot type in the app | Scientific role and concrete use |
|---|---|
| Against diameter | Test/describe size dependence in R2–R6; show whether macrophage or structural relationships vary with object size. Binned trends are descriptive; donor-aware size interactions answer the inferential question. |
| By donor status | Direct ND/Aab-positive/T1D comparisons for R1–R9, displaying zeros and the full spread. Object-level violins/boxes are not independent-donor tests. |
| Donor medians | Show donor support for every hypothesis with one donor per point; pair with donor proportions/count densities where medians would conceal rare events. The card's Kruskal–Wallis test addresses the selected donor summary. |
| Cumulative distribution | R1/R6/R8: separate a change in zero-event frequency from a shift throughout the distribution or an extreme tail. Status-pooled ECDFs weight donors by object number. |
| One curve per donor | R1–R10: show whether a pattern is broadly represented or driven by a donor, including rare macrophage/B-cell profiles. Curves are distributions, not longitudinal paths. |
| Against another measure | R1–R5/R7: examine immune allocation against a second offered object/population measure. The current X selector offers object features and core/peri phenotype proportions, not arbitrary newly derived geometry. Exported joins supply additional comparisons. Object and donor Spearman correlations are distinct from adjusted within-donor effects. |
| Against a donor clinical value | R11 and age context for R8/R9; use donor medians/quartiles and donor-level correlations. Do not replicate the clinical sample size by plotting every object. |
| Density against diameter | R5/R6/R10: reveal the common size–feature combinations behind an overplotted scatter. These are status-normalized frequency contours, not cells/mm² or a separate hypothesis test. |

The **Endocrine objects per donor** card supplies the full Islet/Cluster/Single count and tissue-area-density context for R6. Linked segmentation and single-cell composition panels let the same selected object be inspected across hypotheses. Zeros remain included for the primary view. Axis clipping changes presentation, not the sample; any zero exclusion is a named sensitivity.

### Spatial tab

| Analysis or view | Scientific role |
|---|---|
| Tissue scatter by phenotype | Physical setting for R1–R8 and R10, including extra-islet cells. Tissue position is not UMAP position. |
| Tissue scatter by Leiden | Locate the parent object's composition state in tissue. These colors are not independently discovered cell-neighborhood clusters. |
| Leiden composition | Describe state-defining endocrine proportions; other phenotype bars provide context. Endocrine input features cannot independently validate the clusters. |
| Leiden donor-status and donor membership | Check donor dominance of clusters in R10. These bars show who contributes to a cluster; separately calculate the proportion of each donor's objects occupying each cluster. The two denominators are different. |
| Paired UMAPs with all metadata colorings | Compare composition, immune context, size and donor/clinical associations on the same coordinates. A donor-level color is repeated across that donor's objects. |
| Insulitis prevalence | Operational immune-count threshold prevalence for R1/R6/R8; show the prespecified 15-cell threshold and a limited fixed sensitivity set. The app's aggregate phenotype criterion is not by itself a histopathological diagnosis, and fixed counts depend on object size. |
| Immune cells and beta-cell content | Within-donor paired peri-immune comparison for beta-containing versus beta-poor objects, supporting R1/R8. Keep the chosen beta-proportion cutoff explicit and accompany it with beta presence/counts. |
| Immune gradient around islets | R1/R6–R8: donor-balanced radial enrichment relative to each donor's own outer tissue, with donor curves, counts and support. It separates local concentration from high tissue-wide abundance. |
| Phenotype composition and enrichment bars/heatmap | Core proportions, peri proportions and peri enrichment for all relevant identities; compare cohort and own-tissue baselines. The z value is a descriptive standardized residual, not a hypothesis-test z score. Donor means retain all finite object values, including extreme z scores; the axis-range control changes the bar-chart display only. |
| Distance-to-islet density curves | R1/R6–R8: phenotype or marker-positive spatial distributions with raw-count/proportion modes and both boundary definitions. Status-pooled per-cell curves are descriptive; their height also depends on available tissue and donor cell numbers. |
| Distance-to-islet composition bars | Marker/phenotype positivity or composition among cells in each distance bin; lineage-source interpretation requires the R7 join. Marker rows overlap and are not an exclusive cell-type mixture. ColIV is handled as objects, and Ki67 requires the upstream join. |

### Trajectory and Viewer tabs

| Analysis or view | Scientific role |
|---|---|
| Cell-type-ratio versus pseudotime | Core/peri profiles for R1/R8/R10 on the accepted graph. Export the actual `log1p(cells / max(core endocrine cells, 1))` definition; it is neither a tissue-area density nor a pure cell fraction. Report objects whose denominator required the floor. |
| Marker-expression versus pseudotime | Source-image/scaling inspection only under the current measurement restriction. For biological questions use accepted positive-cell fractions or structure measurements in a derived plot; do not relabel a mean-expression trace as a positive fraction. |
| Raw and rank-normalized DPT | Raw graph ordering and a display with more uniform point spacing; both use the same fitted cohort. Neither measures elapsed time. |
| Overall and donor-status LOESS | Descriptive trends and potential group heterogeneity; use direct donor-aware contrasts for disease differences. Optional plot bands are not automatically donor-clustered intervals. |
| Status/feature pseudotime strip | Describe which group is relatively represented in each bin and the selected feature's bin mean. It does not show a patient's transition between groups. |
| Multi-feature heatmap | Integrate immune composition and external structural measurements for R10. Export bin means alongside row-standardized/clipped display values; row z scores cannot compare absolute abundance across markers. Non-intensity summaries absent from the selector require an exported figure using the same format. |
| PAGA force-directed layout | Show the connectivity of composition profiles and donor contributions. Runtime PAGA Leiden labels are separate from saved composition Leiden labels; apparent branches are not demonstrated biological fates. |
| Linked UMAPs and selected-object drill-down | Follow the same stable object ID from composition space to cells and image for every hypothesis. |
| Whole-slide Viewer and channels | Validate spatial attribution and present representative positive, negative and discordant examples. Display contrast is not a quantitative cross-donor measurement. |

### Implemented modules that are disabled by default

These methods have useful roles if deliberately used on the corrected root. They are not required simply to exhaust a menu, and their presence in code is not evidence that their older phenotype assumptions already support the new data.

| Module and analysis | Role in the revised program |
|---|---|
| Statistics: donor random-intercept mixed model and ICC | Group differences and donor heterogeneity for an appropriate object-level endpoint in R1–R9. Conditional hypotheses require the explicit exposure/covariate specification rather than substituting the module's generic status-only model. |
| Statistics: pairwise contrasts and effect-size forest | Direct disease comparisons with estimates and uncertainty. Keep donor-summary differences distinct from fitted mixed-model contrasts. |
| Statistics: donor-level ANOVA/Kruskal–Wallis and t/Wilcoxon sensitivity | Check dependence on the object-level model using one summary per donor. Record the actual test/selection rule and multiplicity; sparse subsets may be uninformative. |
| Statistics: diameter interaction, fixed-bin tests and Kendall trend | R5/R6 size dependence; distinguish the actual group-by-diameter interaction from within-bin group trends. Ordered ND/Aab/T1D coding is a statistical contrast, not temporal proof. |
| Statistics: trapezoidal AUC over diameter | Descriptive integrated feature-versus-size summary over a common supported diameter range. It is not ROC AUC, classification accuracy, or total biological burden. |
| Statistics: age, sex, autoantibody and adjusted donor models | R8/R9/R11 sensitivity; preserve small donor counts and missingness. Avoid simultaneous adjustment for every recorded attribute. |
| Inflammation: edgeR differential abundance, volcano and table | Donor-count abundance comparison for R1/R8/R9 context, with the actual cell universe, filtering, library normalization and tested phenotype family recorded. Library-normalized counts are relative abundance evidence, not automatically absolute cells/mm²; retain measured tissue-area densities alongside them. |
| Inflammation: distance-weighted infiltration score | Companion to radial and core/peri summaries in R1/R6/R8: weight 1 inside, and 1/max(distance in µm, 1) outside. Keep the numerator's cell universe and chosen denominator explicit for each mode; do not interpret different denominator modes as interchangeable measurements. |
| AI assistant, disabled by default | Navigation or explanation of documented measurements only; it supplies no additional biological evidence or independent statistical test. |

## Comparison and reporting design

1. **One accepted measurement snapshot.** Record the root, code, marker calls, aggregation definitions, object eligibility, tissue areas and representation feature lists. Counts from different repairs cannot be mixed. Use the image/geometry acceptance evidence produced by the ongoing repair; this document does not initiate that repair again.
2. **One coherent endpoint family.** R1–R9 provide the mechanistic association questions; R10 integrates them and R11 is exploratory clinical work. Specify one primary contrast per question before fitting, with marker refinements, regional comparisons and sensitivity analyses in declared secondary families. Multiple plots and repeated representations do not create extra independent confirmations.
3. **Disease comparisons plus the actual conditional test.** Show ND, Aab-positive and T1D summaries; fit the stated exposure relationship and directly test its difference across groups where estimable. Keep pooled, within-group and interaction estimates distinct. Use the already established donor-aware regression, paired comparisons and permutation methods as appropriate to the endpoint; do not force binary, sparse or count endpoints into a generic continuous model just because a tab offers one.
4. **Two complementary scales.** Donor-level summaries establish biological replication; object- and cell-level plots explain heterogeneity. Show counts/densities alongside fractions/enrichment and local measurements alongside whole-tissue context. Unequal donor sampling, spatial dependence and overlapping bands remain relevant after simple aggregation.
5. **Explicit small-population support.** Report donors, objects, positive cells and missing/undefined denominators for each contrast. Selector omission, a zero donor median or a failed model is not evidence of biological absence. Do not replace rare subtypes with a favorable post hoc grouping.
6. **Avoid circular validation.** Record which measurements built each representation. Do not test beta depletion as an independent discovery of a beta-defined cluster or test a fitted macrophage feature as independent validation of DPT. Do not condition a beta outcome indiscriminately on a state defined by that same outcome.
7. **Preserve the alternative outcomes.** A precise null can weaken the specified relationship; a nonsignificant estimate with wide uncertainty cannot. Opposite directions remain reportable. Define any meaningful-effect threshold before looking at results, not after a p-value. This plan does not assign biological branches to the historical results.
8. **Conventional deliverables.** For each question: a direct disease comparison, a donor-support/distribution panel, the relevant paired/conditional relationship, and spatial image examples where the hypothesis is spatial. Use ordinary app-style scatterplots, donor plots, heatmaps, maps and coefficient/interval displays. Keep measured quantities in titles and put cohort, units, definitions and exact statistics in legends. Export the plotted values and stable IDs.

## Availability and derived work

| Availability | Measurements |
|---|---|
| Existing app capabilities, to use on the accepted corrected root | Object/population/structure summaries; donor/clinical plots; endocrine census; core/peri composition; tissue scatter; UMAP/Leiden; DPT/PAGA; enrichment and radial profiles; images and linked drill-down. |
| Existing upstream or research capabilities requiring explicit joins/recalculation on that same root | Ki67 by lineage; marker positivity within a named lineage; structure-qualified cell/overlap fields; the established endocrine adjacency/permutation analysis; existing nerve–vessel measurements; donor clinical models and censored sensitivity. |
| Geometry refinements proposed here, not claimed as app outputs | ColIV boundary gaps and vessel intersections; GAP43-defined nerve-segment neighborhoods; duct-region and B-to-T spatial comparisons; reticular CD35 description; complete region-based small-object niche densities. |
| Dormant code paths requiring compatibility review if used | Statistics, Inflammation and AI. No enablement or deployment change is part of this revision. |

This distinction prevents a missing app selector from discarding usable source data, while avoiding the claim that all desired spatial measurements already exist. Where a measurement is unavailable, its hypothesis remains specified and its status is reported; a different endpoint is not silently substituted.

## Relationship to earlier proposals

| Earlier proposal | Current place |
|---|---|
| Original H1 immune localization and beta quantity | Shared descriptive context; sharper association in R1. |
| H2 matrix and H3 vascular access | R2 and R4. |
| H4 source of LGALS3 | R1 macrophage profiles plus R7 epithelial source. |
| H5 INS content and later composition amendment | Quantitative INS retired; endocrine composition is retained in R5/R10, hormone-negative cells in R5/R9. These are different questions from hormone content. |
| H6 endocrine adjacency | R5, now explicitly evaluated against the fixed-composition spatial reference. |
| H7 neurovascular support | R3 and R4 without an assumed protective function. |
| H8 small endocrine objects | R6 with full-census and tissue-area denominators. |
| H9 disease-stage containment | Direct disease comparison throughout; no longitudinal stage claim. |
| H10 clinical bridge | R11, with any change to the spatial predictor declared in advance. |
| 27 October cycling extension | R9, using the upstream Ki67 measurement and lineage denominator explicitly. |
| 27 October app-organizing principles | R10 plus the full analysis coverage tables above. |

R1–R3 remain the first conceptual priorities because they address myeloid resolution, interface configuration and the additional neural markers most directly. R4–R9 complete the tissue account; R10 and R11 integrate it. The ambition is a connected explanation of spatial heterogeneity, not eleven unrelated claims or a catalogue of significant marker comparisons.

## Evidence and implementation references

The literature comparison is focused, not a systematic proof of novelty. In addition to the sources linked beside the hypotheses, the [27 October proposal](revised_hypotheses_beyond_barlow_20261027.md) records the Barlow panel/method comparison, other relevant pancreas studies and the physiological limits of CD68/CD163/CD4/LGALS3 interpretation. Macrophage–T-cell interactions and early-disease immune changes also appear in [Steenbuck et al.](https://pubmed.ncbi.nlm.nih.gov/42393283/); the refined phenotype/interface relationship, not the general interaction, is the question here.

Current capability and measurement sources inspected for this revision:

- [Plot UI](../app/shiny_app/R/mod_plot_ui.R), [eight plot types and clinical fields](../app/shiny_app/R/plot_types.R), [measure filtering](../app/shiny_app/R/plot_helpers.R).
- [Spatial UI and composition states](../app/shiny_app/R/mod_spatial_ui.R), [neighborhood and distance definitions](../app/shiny_app/R/mod_spatial_server.R).
- [Trajectory UI](../app/shiny_app/R/mod_trajectory_ui.R), [Python DPT/PAGA implementation](../app/shiny_app/python/paga_pseudotime.py), [R trajectory helpers](../app/shiny_app/R/barlow_pseudotime.R).
- [Statistics](../app/shiny_app/R/mod_statistics_server.R), [Inflammation](../app/shiny_app/R/mod_inflammation_server.R), [default feature flags](../app/shiny_app/R/00_globals.R), [Viewer](../app/shiny_app/R/mod_viewer_ui.R).
- [Panel registry](../scripts/uniform/typing/registry/codex-uniform-32__pancreas.csv), [typing rulings](typing_rulings_ledger.md), [historical app-based plan](app_based_hypothesis_plan_20261026.md).

The long user guide and older typing notices contain superseded details. This revision uses the current source for available analyses and the latest owner rulings for biological measurement scope; it makes no claim that the public deployment already serves the local candidate data.
