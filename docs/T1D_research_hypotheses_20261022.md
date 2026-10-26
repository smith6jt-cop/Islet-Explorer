# Ten linked hypotheses for early type 1 diabetes in Islet Explorer

Research proposal dated 22 October 2026. These are candidate hypotheses and analysis specifications, not findings from completed hypothesis tests.

The central question is why endocrine objects exposed to a similar pancreatic inflammatory environment differ in immune penetration, hormone content, and structural support. The proposed program separates tissue exposure, access across an endocrine boundary, and retained endocrine phenotype. It uses the current 20-donor CODEX dataset as a spatial study of early disease. It does not interpret the three donor groups or the app's pseudotime as a longitudinal sequence.

The strongest discovery opportunities are conditional relationships: whether a matrix or vascular configuration changes immune access; whether the cellular source of LGALS3 changes its association with endocrine preservation; and whether small endocrine Clusters occupy different niches from larger objects. Generic claims that T1D has fewer beta cells, more inflammation, or altered small-object abundance already have substantial support, including preliminary findings in the supplied P01 documents.

## Evidence and cohort scope

### The dataset that supports the proposal

The live app selects `data/development/uniform-20donor-20261022-v10/`, built from the uniform processing run under `data/development/uniform-20donor-20261015/`. The latest decisions in [CLAUDE.md](../CLAUDE.md), the actual files, and the [uniform workflow](../scripts/uniform/README.md) govern this proposal. Older 15- and 21-donor research results are not current-cohort evidence. Donor 78 remains excluded.

| Resource | Verified scope | Consequence for inference |
|---|---|---|
| Donors | 20 independent donors: 7 ND, 7 Aab+, 6 T1D | Donors are the biological replicates. Thousands of objects do not create thousands of independent patients. |
| Clinical timing | All six T1D donors have recorded duration 0 to 0.58 years; five have duration 0 | Study recent-onset reserve and injury. Long-duration persistence, progression rate, and remission are outside this cohort. |
| Age | ND 1.58 to 24.60 years; Aab+ 2.20 to 22.00; T1D 3.75 to 25.06 | Age overlaps, but six T1D donors cannot establish age-defined endotypes. |
| Tissue cells | 3,074,596 analyzed cells across 20 images | Tissue-wide spatial context is available beyond the immediate islet ring. |
| Paired analysis objects | 4,763 rows: 1,487 Islets and 3,276 Clusters | Primary core/peri analyses use these paired rows and their eligible cells. |
| Full geometric census | 4,906 Islet/Cluster annotations: 1,491 Islets and 3,415 Clusters | The 143 objects absent from paired tables must be accounted for in density analyses; they are not biological zeroes. |
| Singles | 5,207, counted separately | Under the existing owner rule, use only counts and count densities. Do not give Singles islet-level states, neighborhoods, trajectories, or response models. |
| Embedding | 1,215 objects with at least 10 core cells: 1,202 Islets and 13 Clusters | UMAP/Leiden use four endocrine composition shares. Smaller objects are outside this embedding. |
| Hormones and protein panel | 31 common analyzed channels; all vary | Protein abundance, cell identity, and structural geometry can be studied together. They are not transcript counts or calibrated protein concentrations. |
| Clinical data | C-peptide recorded for 20 donors, with one censored value; HbA1c for 18; Aab flags for 20 | Serum measurements offer an exploratory clinical bridge, with terminal-illness and censoring limitations. |
| Aab+ subgroup | Three single-positive and four multiple-positive donors | Specificity and number are exploratory modifiers. Aab+ alone does not assign clinical stage 1 or 2. |

The analyzed panel is CD31, CD8a, CD3e, SMA, CD4, CD34, HLADR, PDPN, panCK, ECAD, CD163, SST, ColIV, VIM, CD20, LGALS3, B3TUBB, GCG, KRT14, GAP43, CD35, CHGA, PGP9.5, INS, CD44, CD45, NaKATPase, BCatenin, CD68, BActin, and CK19. Ki67 and DAPI are available in source measurements and image channels but are not part of this 31-channel analysis matrix. Other source channels have incomplete donor coverage.

Source compartment-specific measurements undergo REDSEA correction, UniFORM alignment, and a clipped cohort scale before arithmetic aggregation. Preserve each marker's measurement compartment; source Ki67 requires nuclear localization. The assembler, rather than an older `uns` description of log1p measurements, defines this transformation. Current typing and generic scaled-positive flags have different thresholds: `<marker>_pos` uses scaled abundance at least 0.10; endocrine typing uses image-specific hormone levels and a twofold winner margin. A positive INS flag therefore must not substitute for the Beta cell label without an explicit sensitivity analysis.

Several constraints materially affect feasibility:

- There are only 128 typed CD8a T cells in analysis cores; 103 are from two T1D donors. All seven ND donors have zero typed core CD8a T cells. Core CD4 and B cells are rarer still. Primary access tests should use the broad immune union and count-aware models; subtype-specific effects require donor influence checks.
- Exocrine is a residual class, not a validated acinar, ductal, endothelial, or fibroblast label. Core Exocrine and generic Endocrine zeroes arise from the typing cascade. CK19/panCK or VIM/SMA-defined analyses require image-verified operational masks rather than invented cell types.
- Peri bands are tissue-clipped 20 micrometer rings excluding other Islet/Cluster cores. Adjacent bands can overlap. Object-associated counts and counts of unique tissue cells are different estimands.
- Nerve, capillary, and lymphatic annotations are structures measured as areas, not cell phenotypes. The inherited Aab_6505 nerve set was previously rejected; the replacement T1D_6520 nerve set needs fresh adjudication. Lymphatic and capillary annotations also have uneven validation coverage.
- Pancreatic head/body/tail and lobule identities are not validated analysis fields. Fixed spatial tiles can test regional effects, but cannot be called anatomical lobules. Section densities cannot establish whole-organ beta-cell mass.
- The analyzed panel lacks HLA-I/B2M, IAPP, proinsulin, PPY, FOXP3, PD1/TIM3, cytokines, viral proteins, and TCR information. HLADR is class II and does not substitute for HLA-I. Hormone-negative CHGA-positive cells can include unmeasured endocrine identities; they are not identified former beta cells.
- The P01 reports describe Xenium and live-slice function studies, but no usable linked matrices or functional files were found in this project. The Xenium combination notebook contains placeholder paths. Those assays are not inputs to these ten hypotheses.

### How the recent evidence changes the questions

[Barlow et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40232951/) places islets within a coordinated extra-islet immune environment. Its T1D/Aab cases were selected for insulitis, whereas the three ND controls were selected randomly; its states should not be imposed on every donor. [Steenbuck et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42393283/) identifies macrophage–T-cell niches and early endocrine alterations in 88 donors. [Van der Heide et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41673026/) shows strong size dependence and tissue-wide heterogeneity; its analysis excludes objects below 1,000 square micrometers. The local Bodenmiller and Heide PDFs are earlier preprints of these two now-published studies, not independent corroborations.

Recent work also supplies competing interpretations. [Xie et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41477833/) links galectin-3 to autoimmunity through human serum associations and mouse/intervention experiments; it does not establish that all human pancreatic LGALS3 is harmful. [Clarke et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40884792/) and [Dewhurst-Trigg et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41422462/) motivate potentially protective stromal niches, whose precise cell identities cannot be reproduced with this panel. [Richardson et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36696598/) motivates testing neurovascular organization beyond structural density.

The supplied P01 progress report already describes collagen-IV-intact peri-insulitic islets, vessel-adjacent infiltrates, loss of insulin-positive clusters, and functionally abnormal residual beta cells. The added contribution proposed here is to quantify conditional spatial relationships and distinguish their alternatives. It is not a claim to discover these basic observations again.

## The linked aims

| Aim | Hypotheses | Scientific decision |
|---|---|---|
| 1 Explain local immune access | H1 to H4 | Does endocrine injury track a regional immune field, boundary access, structural routes, or source-specific inflammatory niches? |
| 2 Define the pattern of endocrine reserve | H5 to H8 | Is reduced insulin explained by fewer Beta cells, altered abundance within retained cells, endocrine architecture, or selective vulnerability of small Clusters? |
| 3 Connect tissue heterogeneity to clinical phenotype | H9 and H10 | Can immune activity coexist with reserve in Aab+ tissue, and does spatial organization add information beyond endocrine quantity? |

The working dependency is H1 for the spatial exposure definition; H2 and H3 for candidate access constraints; H4 for the cellular inflammatory context; H5 to H8 for distinct reserve readouts; and H9/H10 for donor-level integration. These are analytic dependencies, not a proven biological sequence. A lack of support at one step redirects the linked follow-up tests rather than invalidating the entire program.

## Aim 1 Explain local immune access

### H1 Local immune penetration carries information beyond the pancreatic immune field

**Hypothesis.** Among objects from the same donor and exposed to a similar tissue immune density, greater localization of immune cells inside the core is associated with lower Beta cell density. The association differs from simply having more immune cells in the surrounding tissue.

**Publication link and unresolved question.** Barlow's extra-islet coordination and van der Heide's endocrine abnormalities in immune-poor objects leave open how much local penetration contributes after the regional field is considered. The new contrast is penetration versus exposure, not T1D versus ND immune abundance. [Barlow 2025](https://pubmed.ncbi.nlm.nih.gov/40232951/), [van der Heide 2026](https://pubmed.ncbi.nlm.nih.gov/41673026/).

**Primary test.** Among objects with at least one nearby immune cell, model Beta cell count with core area as an offset, using the within-donor core fraction of the unique core-plus-band immune union as the prespecified exposure. Objects with no nearby immune cells form a separate no-exposure category; their localization fraction is undefined, not zero or successful containment. Report that selection and analyze immune exposure itself in a separate component. Adjust for object size/class, band area, local tissue immune density, and donor effects. Retain separate core and band counts so an apparent ratio effect can be challenged. Use a hurdle/count formulation if zeroes dominate. In a reciprocal localization analysis, condition core immune counts on the nonzero available nearby immune pool. CD3e/CD8a, CD4, CD20/CD35-associated, and CD68/CD163/HLADR strata are secondary, limited by their event counts. Shared band cells require union accounting or explicit overlap correction.

**If supported.** Local access is a candidate discriminator between similarly inflamed objects; field inflammation alone is insufficient. This prioritizes mechanisms addressed by H2 to H4, without demonstrating immune-mediated killing.

**If a meaningful effect is excluded.** Current endocrine depletion may be more strongly related to regional exposure, past injury, or a shared donor process than to current local occupancy. This challenges models that use the amount of visible insulitis as a complete measure of injury.

**Same-data follow-up.** Compare immediate 0–20 micrometer exposure with nonoverlapping 20–100 and 100–250 micrometer tissue annuli, then identify immune-poor/Beta-poor and immune-rich/Beta-rich discordant objects. Refit after leaving out each T1D donor, especially the two dominating core CD8 counts. No current immune cells does not establish that an object was never attacked.

### H2 Matrix continuity is associated with immune exclusion at equal nearby burden

**Hypothesis.** At equal nearby immune burden, gaps in the ColIV-rich endocrine boundary are associated with more core immune localization, while a continuous boundary is associated with peri-islet accumulation.

**Publication link and unresolved question.** The P01 already shows intact ColIV around some peri-insulitic islets. Perivascular/periductal pathology in [Wright et al. 2024](https://pubmed.ncbi.nlm.nih.gov/37881846/) and stromal niches in [Clarke et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40884792/) and [Dewhurst-Trigg et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41422462/) motivate a conditional barrier test. ColIV abundance alone does not measure boundary continuity or tolerance.

**Primary test.** Derive a blinded, image-verified boundary coverage measure from ColIV in a narrow annulus around the endocrine outline, separating nearby capillary basement membrane. Among immune-exposed objects, the primary coefficient is boundary coverage in a model of core immune count with a core-area offset, conditional on peri-islet immune density, perimeter, tissue edge, and vascular area. Handle exposed objects with zero core or zero band counts explicitly; keep objects with no immune cells in a separate exposure component. This avoids placing peri count in both the predictor and an outcome denominator. A coverage-by-peri-density interaction is secondary and asks whether the relationship changes with exposure. VIM, SMA, CD44, ECAD, and BCatenin describe prespecified secondary matrix/adhesion contexts. Exclude structures too small to resolve boundary coverage; report their eligibility separately.

**If supported.** The spatial arrangement of matrix identifies a configuration associated with lower core immune localization beyond total immune exposure. Peri-islet containment becomes a candidate explanation to examine alongside the separate reserve endpoints.

**If a meaningful effect is excluded.** ColIV discontinuity may be a consequence of remodeling, a sectional artifact, or a weak access constraint. This favors vascular topology, field processes, or non-ColIV routes over a simple basement-membrane barrier model.

**Same-data follow-up.** Compare immune-cell distances to observed boundary gaps with positions randomized along the same perimeter; compare gap-centered entry with entry adjacent to CD31-positive structures. Shift the annulus and repeat on raw versus corrected channels. This cannot distinguish proteolytic damage from pre-existing anatomical openings.

### H3 Blood vessel and lymphatic geometry distinguish core and peri immune localization

**Hypothesis.** After conditioning on local vessel abundance and immune density, immune cells localized inside endocrine cores show excess proximity to capillary boundary crossings. Peri-islet immune accumulation has a separable relationship to validated lymphatic structures.

**Publication link and unresolved question.** The P01 reports vessel-adjacent infiltrates. [Medina-Serpas et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41875135/) motivates lymphocyte–endothelial organization, while [Clarke et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40884792/) cautions that vessel-associated niches may also restrain immunity. The test concerns geometry; lymphotoxin signaling and antigen tolerance are not measured.

**Primary test.** Quantify excess proximity of core immune cells to CD31-supported structures crossing the endocrine boundary, relative to spatial nulls preserving tissue masks, object geometry, local immune clustering, and vessel density. Aggregate the observed-minus-expected contact statistic by donor. Use CD34, SMA, and ColIV as complementary structural features. The prespecified alternative analysis tests PDPN-supported lymphatic proximity against peri-islet retention, with blood vessel burden held fixed. Require image review of blood/lymphatic identities; PDPN alone is insufficient.

**If supported.** Candidate anatomical access configurations are associated with immune localization beyond generic vascular enrichment. This nominates a structural explanation to investigate alongside H2.

**If a meaningful effect is excluded.** Vessel abundance may reflect endocrine support, tissue contraction, or general anatomy rather than preferential immune entry. A lymphatic-only association would redirect attention to retention or drainage-associated organization, not prove either process.

**Same-data follow-up.** Separate capillary–ColIV gaps, vessel-adjacent HLADR-positive immune-negative candidate stroma, and peri-lymphatic B/T aggregates. Repeat with equal-area reference tissue and with the alternative structural masks already retained in the source run. Static 2D images cannot demonstrate the direction of trafficking, perfusion, or lymph flow.

### H4 The cellular source of LGALS3 determines its association with endocrine injury

**Hypothesis.** At equal combined myeloid/epithelial LGALS3 signal and immune burden, a greater myeloid allocation of that signal is associated with reduced Beta cell density. Endocrine-source LGALS3 and immune penetration are separate secondary readouts.

**Publication link and unresolved question.** [Xie et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41477833/) provides a recent galectin-3 mechanism in experimental autoimmune diabetes, whereas [Steenbuck et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42393283/) places myeloid–T-cell neighborhoods at the human islet interface. The unresolved human tissue question is source and context, not whether galectin-3 is uniformly pathogenic.

**Primary test.** Derive LGALS3 within typed macrophages and image-verified CK19/panCK/ECAD-positive epithelial masks. Define the source allocation as myeloid LGALS3 divided by the combined myeloid-plus-epithelial LGALS3 signal. Here signal means the sum of scaled cell means, not integrated protein mass. Test its within-donor association with Beta cell density, conditional on the combined signal per tissue area, source-cell densities, and local immune density. Treat a near-zero signal pool as unassessable, not an arbitrary ratio; fix the minimum pool from background reliability. Analyze endocrine-source signal separately. CD68, CD163, HLADR, CD45, CD44, and VIM are secondary descriptors. Use cell-area/overlap and neighborhood checks to distinguish intracellular source assignment from adjacent staining.

**If supported.** Whole-islet or whole-tissue LGALS3 averages obscure biologically different niches. Myeloid-source localization would nominate a more specific pathway for subsequent mechanistic work.

**If a meaningful source contrast is excluded.** LGALS3 may primarily report the number of source cells, general epithelial remodeling, or a shared response to injury. A reversed association would challenge transferring a harmful myeloid interpretation to all tissue compartments.

**Same-data follow-up.** Match high-LGALS3 neighborhoods on cell composition; compare macrophage-rich/epithelial-rich regions and core/peri/distal tissue; examine CD163 and HLADR jointly without assigning binary M1/M2 function. Treg suppression, secretion, and galectin-3 treatment response are not established by these markers.

## Aim 2 Define the pattern of endocrine reserve

### H5 Reduced insulin signal reflects both cell quantity and the phenotype of retained cells

**Hypothesis.** Within donors, reduced insulin abundance in inflamed niches is not fully explained by fewer typed Beta cells; retained Beta cells have lower insulin abundance, with a distinct contribution from CHGA-positive hormone-negative cells.

**Publication link and unresolved question.** [Steenbuck et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42393283/) and [van der Heide et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41673026/) separate endocrine molecular changes from complete cell depletion. This panel can test insulin-content decomposition, but cannot reproduce their IAPP/proinsulin endpoints or prove the functional defect already described in the P01.

**Primary test.** For a prespecified reviewed core-cell set, use the additive identity `sum(INS_all) = N_Beta * mean(INS_Beta) + sum(INS_nonBeta)`. INS denotes a scaled per-cell mean; this sum is not integrated fluorescence or protein mass. Retain the non-Beta contribution, including measured unresolved cells when the source eligibility definition admits them. Test the association between local immune penetration and within-Beta mean INS, conditioning on Beta cell density, object size, and image background. A median is a robust secondary endpoint, not an additive decomposition; an object with no Beta cells has an unavailable within-Beta mean. Quantify CHGA-positive hormone-negative cells separately. NaKATPase, BActin, ECAD, BCatenin, CD44, and lineage-conditioned HLADR are prespecified secondary state/QC readouts. Do not divide by BActin or NaKATPase without demonstrating that they are stable references.

**If supported.** Quantity alone misses a spatially structured hormone phenotype in retained cells. This yields a candidate tissue phenotype for future functional linkage.

**If a meaningful within-cell effect is excluded.** Cell depletion and composition may account for the measured difference among selected typed Beta cells. The result does not exclude insulin attenuation in cells that failed the typing gate. The source residual and independent-mask analyses are needed before extending the conclusion to all endocrine cells.

**Same-data follow-up.** Because INS contributes to typing, repeat across hormone thresholds and winner margins; inspect continuous INS within an independently reviewed endocrine mask; compare CHGA-supported cells, unresolved-margin cells retained in source audits, and nearby alpha/delta cells. Repeat with raw and corrected measurements and exclude weak staining acquisitions in sensitivity analyses. A robust hormone-negative population still cannot establish beta lineage, dedifferentiation, or transdifferentiation, particularly without PPY.

### H6 Alpha and delta topology modifies local Beta cell preservation beyond composition

**Hypothesis.** At the same endocrine composition and immune exposure, objects with greater alpha/delta adjacency to Beta cells retain more Beta cells or higher within-Beta INS abundance than expected from random placement of their existing cell types.

**Publication link and unresolved question.** Size-dependent endocrine architecture in [van der Heide et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41673026/) makes raw hormone proportions an inadequate test of spatial support. The new question is whether alpha/delta contact organization adds information after their abundance is held fixed.

**Primary test.** Within sufficiently large eligible objects, measure excess Beta–Alpha and Beta–Delta adjacency over label permutations preserving counts and geometry. Use one prespecified combined contact statistic and test its adjusted main association with within-Beta INS abundance at fixed immune exposure; Beta density and contact-by-immune interaction are secondary endpoints. Compare matched composition/size strata. Sparse delta counts should widen uncertainty rather than justify unstable object-specific estimates. Restrict the primary test to objects with enough cells of the required types.

**If supported.** Islets with identical composition can have different reserve-associated architecture. This identifies a spatial modifier for local endocrine vulnerability.

**If a meaningful effect is excluded.** Composition or shared tissue conditions may explain the apparent support relationship; adjacency adds little at the resolution measured here. A reversed effect would suggest remodeling around residual Beta cells rather than simple buffering.

**Same-data follow-up.** Repeat using true cell-boundary contacts versus centroid neighborhoods, separate alpha and delta contributions, and condition on vessel/nerve proximity from H3/H7. The analysis does not demonstrate glucagon or somatostatin signaling, functional protection, or direction of causality.

### H7 Neurovascular alignment identifies supported endocrine niches beyond density

**Hypothesis.** At equal nerve and capillary abundance, excess nerve–vessel alignment near endocrine objects is associated with greater endocrine reserve or higher retained-cell hormone abundance.

**Publication link and unresolved question.** [Richardson et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36696598/) and [Campbell-Thompson et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33753784/) show why structural preservation and function should be distinguished. Barlow's extra-islet structural findings motivate testing organization beyond the density changes already reported. [Barlow et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40232951/).

**Primary test.** First derive a common, INS-independent B3TUBB fiber mask using morphology and blinded review, with GAP43 as corroboration where present. This is essential because the original 15 donors' inherited nerve classifier uses INS as well as B3TUBB, whereas additions use B3TUBB thresholds; using those masks directly against an INS outcome risks circularity. Quantify fiber proximity/alignment to CD31-supported structures in the peri region. Test the donor-specific association of excess-over-null alignment with within-Beta INS abundance, adjusting for separate nerve/vessel abundance, area, perimeter, and immune exposure. Nulls must preserve tissue support and local structural abundance. PGP9.5 requires exclusion of endocrine staining; inherited and replacement masks are secondary comparisons.

**If supported.** A reserve-associated spatial configuration distinguishes endocrine niches despite similar structural density, nominating a support mechanism for further investigation.

**If a meaningful effect is excluded.** Neurovascular alignment adds little independent hormone-phenotype information at the measured precision. Separate density/geometry follow-ups can then test whether apparent structural preservation is compatible with tissue contraction or passive scaffold persistence; an alignment null alone cannot establish those alternatives.

**Same-data follow-up.** Compare core/peri/distal fiber area, length, and alignment; challenge each result with B3TUBB versus GAP43-supported masks and omission of uncertain images. If density effects vanish after geometry adjustment but alignment remains, prioritize organization. This panel cannot establish sympathetic identity, neurotransmission, reinnervation, or glucagon counterregulation.

### H8 Small Beta bearing Clusters are selectively scarce in particular tissue niches

**Hypothesis.** The scarcity of small Beta-bearing Clusters in T1D is concentrated in immune-rich, epithelial/stromal remodeling fields rather than distributed evenly across available tissue contexts. Define small by 2–9 typed endocrine cells for the primary analysis so immune entry itself cannot move an object into a larger size category.

**Publication link and unresolved question.** [Murrall et al. 2025](https://pubmed.ncbi.nlm.nih.gov/41223263/) reports reduced small extra-islet insulin-positive objects, while the whole-organ study of [Lehrstrand et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42172326/) cautions against universal absence across clinical settings. The latter is a single late-onset pancreas and cannot be treated as a direct comparator to this recent-onset cohort. The new question extends the local P01's cluster-loss observation to a specific niche-by-size interaction.

**Primary test.** Count Beta-bearing Clusters with 2–9 analyzed, typed endocrine cells per tissue area in a prespecified fixed spatial grid, classifying fields using immune burden and image-verified epithelial/stromal markers independently of the Beta endpoint. The primary contrast is the disease-by-field association on this small-Cluster density. Comparisons with small and larger Islets are secondary and cannot isolate size from object class. Only 13 paired Beta-bearing Clusters have at least 10 typed endocrine cells, including one T1D object, so a within-Cluster size interaction is not a supported primary test. Verify these counts against the full source census before fitting. Repeat the small-object definition using reviewed geometric area and raw whole-cell counts because endocrine depletion and hormone gating can also change endocrine-count size. Account for all 143 geometric objects missing from paired rows. Objects with entirely excluded/unmeasured cores have uncertain Beta identity, not established Beta negativity; include their geometry and show classification/exclusion bounds. Analyze the full eligible source census, not just embedding members. Standardize grid origin/scale and tissue-edge handling before inspecting outcome maps.

**If supported.** Small-object scarcity is anatomically selective and points to a local vulnerability or maintenance niche beyond generic small-islet loss.

**If a meaningful interaction is excluded.** More spatially uniform scarcity, developmental abundance, section sampling, or detection limits may explain the distribution better than the proposed niche. This would constrain a niche-selective mechanism without refuting small-object depletion itself.

**Same-data follow-up.** Compare 100, 200, and 400 micrometer grids; hormone composition and CHGA-positive hormone-negative counts in eligible Clusters; alternative size/eligibility definitions; and age-matched donor subsets. Use Singles only as donor-level counts/densities alongside this analysis. Source Ki67 can provide a separately validated exploratory proliferation-associated readout after nuclear localization and background review; positivity cannot prove regeneration or lineage conversion.

## Aim 3 Connect tissue heterogeneity to clinical phenotype

### H9 Aab positive tissue can show immune exposure with retained local reserve

**Hypothesis.** At comparable regional immune burden and donor-wide Beta cell quantity, matrix continuity is more strongly associated with low local immune penetration in Aab+ than recent-onset T1D tissue. This tests a spatial containment configuration beyond the expected group difference in overall Beta abundance.

**Publication link and unresolved question.** The broad single-/multiple-Aab design of [Steenbuck et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42393283/) and the architecture in [van der Heide et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41673026/) motivate separating exposure from injury. [Torabi et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37246834/) supports examining immune heterogeneity without assuming age alone defines an endotype.

**Primary test.** Estimate the Aab+-versus-T1D difference in the H2 boundary-coverage association with core immune density; this group-by-coverage contrast is the single primary test. Compare donor-specific within-donor coverage slopes using a small donor-level model accounting for donor-wide Beta quantity, or an equivalent hierarchical slope model. A donor-fixed intercept absorbs donor-constant main effects; it does not identify a separate coefficient for donor-wide Beta quantity. Require joint overlap in Beta quantity, boundary coverage, and immune exposure, along with sufficient event-bearing donors. If these conditions fail, the adjusted group contrast is unestimable rather than negative; report descriptive within-donor relationships without extrapolating. Use the exposure/localization separation in H1/H2: no nearby immune cells is not successful containment. H3/H7 support features are secondary alternatives, not interchangeable primary exposures. Describe fields combining local reserve with low penetration using cutoffs fixed by measurement reliability and training-donor ND reference data. Keep every tissue tile in the density census, including object-free tiles; penetration/boundary features are structurally unavailable there, not evidence of successful containment. Report these negative spaces separately. ND supplies a reference distribution, with age/size sensitivity analyses. Do not force the three groups onto an ordered score.

**If supported.** Matrix continuity has a different association with local immune exclusion across clinical groups beyond the supported range of donor-wide Beta quantity. Test the relationship to reserve-associated exposed fields in the secondary analysis; this primary result alone does not establish endocrine sparing or predict an individual's progression.

**If a meaningful difference is excluded.** The matrix-associated containment relationship may be shared between groups or weak in both. Distinguish those outcomes before asking whether total exposure, endocrine quantity, stage ascertainment, age, or acquisition better explains the group distinction.

**Same-data follow-up.** Examine the three single- versus four multiple-Aab donors descriptively, then antibody specificity and broad B/T/myeloid neighborhood composition. Review CD20-positive/CD45-positive candidate B cells, CD35 localization, and HLADR source on images. Probe age continuously in all groups rather than discovering age thresholds in six T1D donors. No TLS, antigen specificity, Treg state, or clinical stage is inferred from these limited markers.

### H10 Spatial niche information improves the clinical interpretation of endocrine quantity

**Hypothesis.** A compact score of immune access and structural support adds information about donor C-peptide beyond sampled Beta cell quantity, age, and glycemic context alone.

**Publication link and unresolved question.** Human studies distinguish endocrine quantity, altered phenotype, and spatial context, yet their connection to retained secretion remains incomplete. [Steenbuck et al. 2026](https://pubmed.ncbi.nlm.nih.gov/42393283/), [Dewhurst-Trigg et al. 2026](https://pubmed.ncbi.nlm.nih.gov/41422462/). Fixed-tissue features can be tested as correlates of recorded C-peptide; they cannot establish secretion by an imaged object or treatment responsiveness.

**Primary test.** Compare donor-held-out predictive log scores for two small, prespecified models of log C-peptide. Uncensored observations contribute their predictive log density; the donor recorded as below 0.02 ng/mL contributes the log predictive probability below that bound, not squared error to 0.02. Fit a strongly regularized predictive distribution and its variance only on training donors. The baseline uses sampled Beta cell quantity, age, and one glycemic context measure. The augmented model adds one direction-specified access/support score from H1–H8, with inclusion fixed by blinded measurement reliability rather than whole-cohort association results. Any learned selection, weights, scaling, imputation, or tuning must be repeated entirely within training donors. Report the augmented-minus-baseline score difference with donor-level uncertainty. Use a baseline-conditioned parametric benchmark preserving covariates and documented censoring rules for incremental information; a raw outcome permutation would test global association, not this conditional claim. If assay censoring rules are insufficiently documented, present that benchmark as a sensitivity analysis rather than a calibrated significance test. Avoid using disease group as a shortcut, and report a disease-group-conditioned sensitivity to reveal whether the score mainly distinguishes groups.

**If supported.** Spatial organization may improve interpretation of residual endocrine quantity and nominate a tissue phenotype for future prospective validation. The current result would remain an exploratory clinical association in 20 donors.

**If a meaningful improvement is excluded.** Composition/quantity may be sufficient for this cohort's recorded C-peptide, or terminal physiology may limit the value of this clinical endpoint. A complex spatial score would not be justified for clinical interpretation here.

**Same-data follow-up.** Separate failure of the clinical bridge from failure of the tissue mechanism: test the fixed score against H5's within-cell INS phenotype in held-out donors; examine HbA1c, admission/peak glucose, BMI, and procurement variables as sensitivity readouts. Report group-specific held-out scores, show all six T1D observations individually, and treat a T1D-only model as descriptive. C-peptide at organ donation is not a standardized stimulated test and may be affected by glucose, treatment, renal function, and terminal illness; the available metadata do not remove all such influences.

## Analysis rules that make either outcome useful

### Biological replication and model size

Use within-donor contrasts for the mechanistic hypotheses, followed by donor-level estimates and uncertainty. Decompose object-level exposures into their donor mean and within-donor deviation so between-donor acquisition differences do not masquerade as local biology. Keep group contrasts at the donor level. Donor and disease are partly inseparable from acquisition; neither a random intercept nor Harmony proves those effects have been removed.

Use a minimal prespecified covariate set. Size, perimeter, and tissue/band area belong where the geometry requires them; age and sex enter donor-level sensitivity analyses when support permits. Do not simultaneously fit every clinical/handling field with 20 donors. Report age-overlap subsets, leave-one-donor-out stability, and original-acquisition/addition sensitivity where possible. ICU time, recovery type, transit time, cause-of-death categories, RIN, and pancreas weight are handling/confounding checks, not a fishing expedition for predictors. Preserve missing values and the C-peptide censoring bound; never turn missing HbA1c into zero.

Sparse outcomes require count-aware or hurdle models, not an assumption of Gaussian abundance. With separation or few event-bearing donors, show event counts and donor-level intervals rather than unstable interaction estimates. Donor bootstrap/permutation calculations must resample or permute donors and repeat feature fitting as required. A very large number of object-level observations cannot rescue a contrast supported by only one or two donors.

### Definitions and spatial controls

Keep cell count, cell density, phenotype share, marker abundance, structure area, boundary coverage, and excess proximity as distinct quantities. Use denominators appropriate to the question; loss of Beta cells can inflate every other endocrine share. Retain absolute density and tissue-area companions to compositional results. A panel-wide positive flag is not a universal biological threshold.

Use the current 20 micrometer peri band as the shared reference, then derive 10/30 micrometer and distal-region sensitivity analyses offline from existing tissue/cell geometry. Shared bands need deduplication for union endpoints or overlap-aware analysis. Existing `min_dist_*` summaries are distances to mean core-cell coordinates, not boundary/contact distances; existing enrichment z-scores use pooled cohort frequencies rather than donor-local nulls. Recompute the actual distances and null models specified here. Link source tables by persistent donor/object_id, never row order or concatenated duplicate peri records. Treat missing/invalid structure masks as unassessable, never as zero structure. Reject or review uncertain masks without reference to the desired hypothesis direction. Analyze tissue edges, object fragmentation, and cell spillover explicitly.

The app's threshold of 15 broad immune cells in core plus band is an operational display. It must not be equated automatically with the full histopathological consensus definition of insulitis. Use continuous counts/densities and prespecified threshold sensitivities for the primary research questions.

### Circularity and interpretability

The active UMAP/Leiden representation is built from Beta, Alpha, Delta, and CHGA-positive hormone-negative shares. Its separation cannot independently validate those same shares. The app's Barlow-style diffusion pseudotime uses core/peri composition, an ND root, and an orientation chosen so Beta ratio decreases. An association of Beta ratio with that pseudotime is partly built in. Use embeddings for inspection and example selection, not temporal or mechanistic validation.

INS-based typing creates selection/truncation when testing within-Beta INS. Challenge that dependence using independent reviewed endocrine masks, threshold variation, and continuous measurements. For any trajectory-based secondary test, remove the endpoint and its defining/gating features from the fitted representation; fit only within training donors when predicting held-out outcomes. Keep a record of phenotype, segmentation, and feature provenance.

### Support rejection and uncertainty

For each primary contrast, set a smallest biologically meaningful effect before fitting the outcome model. Choose it using assay reproducibility, plausible tissue consequences, and blinded image review; do not declare a universal threshold across count, contact, and protein outcomes. Simulate precision under the observed donor/event counts before interpreting a negative result.

Classify results as support, evidence against the prespecified meaningful effect, or inconclusive. A nonsignificant p-value is not evidence that the hypothesis is false. A meaningful negative result requires an interval narrow enough to exclude the relevant effect, or a reproducible opposite association. With six T1D donors, some hypotheses will remain inconclusive even if technically measurable. The value of the specified negative branches is that they discriminate alternatives when precision is adequate; importance cannot be guaranteed by statistical nonsignificance.

Treat the ten primary tests as one declared multiplicity family if they are promoted to a formal testing program. Use BH control across all ten slots, with an explicit untestable status for unavailable endpoints and no fabricated p-values. Apply separate declared families to secondary marker modules. These proposals follow inspection of the app and data inventory; they are not preregistered independent confirmations. Freeze the analysis plan before running the tests and describe the results accordingly.

## Coverage of all data components

| Data component or marker module | Main role | Hypotheses |
|---|---|---|
| INS, GCG, SST, CHGA | Endocrine identity, quantity, hormone phenotype, architecture | H1, H5, H6, H8, H9, H10 |
| CD45, CD3e, CD8a, CD4 | Broad immune access and sparse T-cell secondary analyses | H1, H2, H3, H4, H9 |
| CD20, CD35 | Image-verified B-associated/complement-receptor context; no automatic dendritic or TLS identity | H1, H3, H9 |
| CD68, CD163, HLADR | Myeloid context and source-specific class-II-associated signal | H1, H3, H4, H5, H9 |
| LGALS3, CD44 | Source-conditioned inflammatory/remodeling state | H2, H4, H5, H8 |
| ColIV, SMA, VIM | Boundary organization and operational stromal context | H2, H3, H4, H7, H8, H9 |
| CD31, CD34, PDPN | Reviewed blood/lymphatic structure geometry | H2, H3, H7, H9 |
| B3TUBB, GAP43, PGP9.5 | Fiber morphology and neurovascular alignment with endocrine-staining controls | H7 and H6 follow-up |
| CK19, panCK, KRT14, ECAD, BCatenin | Image-verified epithelial fields and adhesion/remodeling descriptors | H2, H4, H5, H8 |
| NaKATPase, BActin | Retained-cell membrane/structural phenotype and stability checks | H5; H10 secondary modules |
| Source Ki67 and DAPI | Nuclear localization, assay QC, exploratory proliferation-associated follow-up | H8; quality checks across aims |
| Full tissue, object polygons, cell coordinates, signed distances | Core/peri/distal fields, contacts, tissue-edge/null models | H1 through H9 |
| Cell/nuclear morphology and 12 retained mean morphology features | Segmentation, contact, size, and staining robustness checks | H2, H5, H6, H7, H8 |
| Census and object classes, including permitted Single counts | Density, selection accounting, small-object context | H8; size adjustment across aims |
| Embeddings and 31-channel heatmaps | Exploratory inspection, donor dominance, discordant-example review | All aims; no independent validation of input features |
| C-peptide, HbA1c, glucose, BMI, Aab specificity/count, age/sex | Clinical bridge and heterogeneity | H9, H10; sensitivity analyses |
| Procurement and tissue metadata | Handling/sampling alternatives | All donor-level comparisons |

Coverage means each component has a defined biological or quality-control purpose. It does not require including every marker in each model or interpreting every available metadata field as a mechanism.

## Priority and concrete deliverables

Start with H1, H4, H5, and H8: they use the project's distinctive tissue context, cellular LGALS3 source, endocrine phenotype, and small Clusters. H2/H3 follow after boundary/vascular review; H7 requires nerve-mask adjudication. H6 needs sufficient cell-type counts and contact reliability. H9 is an exploratory seven-versus-six donor comparison. H10 is a final exploratory integration after the contributing modules have stable definitions; do not optimize its score to make it positive.

Each hypothesis should produce a donor-level effect plot, an object/tile-level descriptive map, eligibility/event counts, representative image overlays chosen under a fixed rule, uncertainty and multiplicity results, and the specified alternative test. Publish all prespecified outcomes, including untestable or inconclusive ones. Retain the full tissue/annotation audit beside the results so apparent biological absence can be distinguished from missing or rejected measurements.

### Project source locations

- [Current cohort and processing decisions](../CLAUDE.md) and [uniform workflow](../scripts/uniform/README.md).
- [Source grant](../literature/PO1-Immune%20Function%20and%20the%20Progression%20to%20Type%201%20Diabetes.docx), Project 1 aims; [progress report](../literature/P01%20Project%201%20RPPR.docx), Aim 1/Figure 2, Aim 2, and Aim 3. These provide context and known preliminary findings, not additional served assays.
- Local literature: [Barlow](../literature/Barlow2025.pdf), [Steenbuck preprint](../literature/Bodenmiller2025.pdf), [van der Heide preprint](../literature/Heide.pdf), [Richardson](../literature/Richardson2023.pdf), and [Campbell-Thompson](../literature/Campbell-Thompson2021.pdf). Use the linked 2026 publications above when citing the final Steenbuck and van der Heide studies.
- Current source cells: `data/development/uniform-20donor-20261015/cells/`; aligned continuous measurements: `typing/aligned/` below that run; typing calls/audits and structural objects are retained beside them.
- Current served summaries, cell flags, signed distances, raw images, and drill-down files: `data/development/uniform-20donor-20261022-v10/`.

The literature review is targeted to the ten hypotheses, not a systematic review establishing worldwide novelty. Donor overlap with external nPOD studies must be checked before claiming independent biological replication. All proposed follow-ups above use this project's existing images, measurements, geometry, or clinical tables; none requires assuming an unprovided transcriptome, lineage assay, or functional experiment.
