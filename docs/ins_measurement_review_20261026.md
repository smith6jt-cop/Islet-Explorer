# INS measurement review and endpoint amendment

The owner states that the CODEX images were not acquired for quantitative comparison of INS signal. That acquisition statement governs this dataset. This review concerns endpoint validity and plan changes; it does not interpret the hypothesis-test results.

## Published evidence

**A relevant CODEX example exists.** Iwasaki et al., *JCI Insight* 2026, “p21-senescent cells drive pancreatic islet dysfunction through targetable paracrine signaling in type 2 diabetes,” report INSULIN and C-PEPTIDE protein comparisons among p21-defined beta-cell groups in human pancreas. Their spatial analysis combines CODEX and Xenium. It is a precedent for protein-signal analysis on the platform, not a demonstrated match to the acquisition conditions of this dataset. [PubMed](https://pubmed.ncbi.nlm.nih.gov/42262869/) · [Article results](https://insight.jci.org/articles/view/197310/version/2).

The supplement describes a 38-marker PhenoCycler-Fusion protocol and QPTIFF images. The verified methods do not establish matched exposure, saturation behavior or a comparable cross-image INS measurement scale for our images. Consequently this reference does not resolve the owner's stated acquisition limitation. [Supplementary methods, PDF pages 18–19](https://insight.jci.org/articles/view/197310/sd/pdf/render/1).

**Barlow provides a count-based precedent.** The published CODEX study reports antibody/exposure optimization, INS exposure of 200 ms and proinsulin exposure of 40 ms. Marker values are z-normalized within donor for classification, and beta cells are identified using proinsulin. The relevant endocrine comparisons use beta-cell counts/frequencies and presence/absence categories; they do not establish a retained-beta-cell INS-abundance endpoint for this dataset. [Barlow et al., eLife 2025](https://pubmed.ncbi.nlm.nih.gov/40232951/); local methods in `literature/Barlow2025.pdf`, pages 19–21.

**The previously cited quantitative hormone study uses another platform.** Steenbuck et al. use imaging mass cytometry. Their INS/proinsulin/IAPP analyses cannot by themselves establish validity of quantitative fluorescence measurements from this CODEX acquisition. [Published study](https://pubmed.ncbi.nlm.nih.gov/42393283/). The supplied `~/Documents/papers/Bodenmiller2025.pdf` is its earlier manuscript, with acquisition and analysis methods on pages 25–28.

The local human-pancreas atlas also uses CODEX protein measurements with normalization, but the inspected endocrine comparisons use CPEP/IAPP/SST and other proteins; its INS-high transcript populations must not be confused with a CODEX INS endpoint. The inspected Heide manuscript uses MICSSS immunohistochemistry. Neither provides an acquisition-matched justification for the current INS endpoints.

## Revised endpoint scope

These are methodological amendments. No new model fits have been run in this review. Existing intensity-based computations remain historical, with no new biological interpretation assigned.

| Analysis | Amendment | What this measures |
|---|---|---|
| H5 | Remove mean/median INS, summed-signal decomposition and raw/corrected intensity variants from the current hypothesis endpoints. Retain existing beta-cell and CHGA-positive hormone-negative counts/densities as population measurements. | Cell quantity and categorical composition. These do not test insulin abundance within retained beta cells. A new composition hypothesis must be stated explicitly rather than relabeling these counts as the original H5 test. |
| H6 | Use the already specified beta-cell-density endpoint and retain fixed-composition adjacency/permutation measurements. Keep the existing count-endpoint analysis separately identified from the historical INS models. | Association between measured endocrine architecture and label-defined beta-cell density. The current graph eligibility requires at least two beta cells, so beta-cell presence cannot be an outcome in this selected graph set. |
| H7 | Replace the INS response with beta-cell presence and beta-cell count/density using the existing presence/count approaches. Retain observed nerve–vessel proximity and separate structure/geometry covariates. | Label-defined beta-cell quantity/presence in relation to proximity. Include eligible objects without beta cells: 2,167 objects have measured proximity, of which 1,221 lack typed beta cells; the previous INS analysis used 946 beta-containing objects. No new image extraction is required for these endpoints. |
| H10 | Retain clinical C-peptide models. Remove the tissue-INS prediction follow-up from current use. | The clinical outcome is an independently recorded measurement. Replacing tissue INS with the exact beta-density predictor already in the baseline would duplicate a predictor and is not an equivalent follow-up. |
| H1/H2/H3/H8/H9 | No direct continuous INS endpoint requires replacement. | Existing label/count/geometry endpoints remain separate from quantitative INS abundance. |
| H4 | INS enters through accepted beta-cell labels, not the continuous INS response. | Continuous LGALS3 allocation is a separate measurement assumption. The INS-specific acquisition statement is not automatically extended to LGALS3 or other markers. |

The original H1–H10 plan remains a historical document. Comparisons across disease states remain the primary presentation. Existing accepted cell identities, counting corrections, edge exclusions, donor exclusions and graph/image measurements are reused.

Raw/corrected variants, z-scoring, rank transforms or within-image comparisons will not be treated as automatic solutions to the stated acquisition limitation. This amendment does not introduce a new calibration project or manual-review prerequisite.

## Existing export

The previously delivered `output/prism_export_20261026/` package contains historical INS-based analyses. Its measurement notice identifies the affected endpoints. Those plotted numbers are preserved for provenance; their presence in the package does not make them current quantitative INS tests. New endpoint results have not been substituted into those tables.
