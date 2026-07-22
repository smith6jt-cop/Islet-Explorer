# CODEX first-pass RESTORE pairs — empirical checks

Measured on `results/old_data_plots/all_data.tsv` (2,676,466 cells, all 15 donors, raw QuPath
`Cell: <marker>: Mean` in native 16-bit units).

**Caveat on method.** "Counterpart-high" below is the top decile of the counterpart within the
donors that carry it. That is a *proxy* for RESTORE's negative control, not RESTORE itself —
RESTORE fits two components on the co-expression cloud and takes the cluster with maximum spread
on the reference axis. These numbers rank candidate counterparts; they do not predict the divisor.

---

## 1. The `>50` floor does not separate background from signal

Requested check: does 50 sit between background and signal in native units? **No.** It sits
inside the background for every abundant marker.

| marker | median | % of all cells > 50 |
|---|---:|---:|
| DAPI | 3999 | 99.99% |
| CHGA | 290 | 99.58% |
| ECAD | 7633 | 99.24% |
| panCK | 587 | 96.10% |
| CD45 | 132 | 80.89% |
| GAP43 | 78 | 68.00% |
| CD31 | 63 | 60.09% |
| VIM | 73 | 59.91% |
| PGP9.5 | 54 | 53.26% |
| B3TUBB | 45 | 44.15% |
| SMA | 30 | 28.05% |

**Consequence:** the floor is effectively a zero-drop for the five gate markers, not a background
floor. It only bites the dim markers (SMA, B3TUBB) — which are not first-pass gates. It is
harmless here, but it is not doing the job its name implies, so do not rely on it.

---

## 2. `epCAM` is not viable cohort-wide → substituted with `panCK`

`epCAM` is acquired on **8 of 15 donors** (6356, 6450, 6479, 6505, 6516, 6538, 6550, 6563);
non-null on 53.1% of cells. The 33-channel donors (112, 6521, 6533, 6534, 6548, 6549, 6551) lack
it. A counterpart present on half the cohort cannot produce a comparable cohort-wide positivity
call, and per-subset stratification would mean two different background definitions for PGP9.5.

Candidate counterparts, scored on the two properties that matter — PGP9.5-cleanliness of the
negative control, and endocrine contamination. Restricted to the 8 donors carrying epCAM so the
comparison is like-for-like (n = 142,171 per set):

| counterpart | PGP9.5 median | PGP9.5 p90 | % PGP9.5 > 500 | % of set that is endocrine-high |
|---|---:|---:|---:|---:|
| ECAD | 9.3 | 50.0 | **0.87%** | 20.26% |
| panCK | 14.4 | 62.8 | 1.46% | **18.14%** |
| CK19 | 15.8 | 67.1 | 1.73% | 17.48% |
| epCAM | 14.9 | 80.3 | 4.06% | 22.80% |
| *(all cells)* | 36.2 | 174.3 | 6.22% | 10.00% |

**`panCK` is chosen** — exocrine-restricted like epCAM, present on all 15 donors, and on both
metrics a *cleaner* negative than epCAM (2.8× less PGP9.5-bright contamination, less endocrine).

Cohort-wide on all 15 donors, panCK-high cells (n = 267,649) have PGP9.5 median 29.8 vs 53.5 for
all cells. Per-donor median within panCK-high spans 4.3–94.9 (**21.9×**) — large donor variation,
which is what per-image RESTORE thresholds exist to absorb, but it means a pooled threshold would
be badly wrong.

### Two stated premises the data does not support

1. **"epCAM-high isolates exocrine, which is definitively PGP9.5−"** — epCAM-high is the *most*
   PGP9.5-contaminated of the four candidates (4.06% > 500, vs 0.87% for ECAD) and the *most*
   endocrine-enriched (22.80%). The substitution to panCK is therefore an improvement on the
   original intent, not a compromise forced by availability.

2. **"ECAD includes endocrine, which is PGP9.5+, so it would seed the background with true
   neural signal"** — ECAD-high *is* endocrine-enriched (20.26% vs 10% baseline), so the first
   half holds. But its PGP9.5 distribution is the **cleanest of all four** (median 9.3, 0.87%
   above 500). The endocrine cells captured in ECAD-high are not PGP9.5-bright in this data at
   this threshold. So the mechanism the rule guards against is not visible here.

   **ECAD is nonetheless not recommended for PGP9.5**, for an independent reason: ECAD is already
   the counterpart for CD45 and VIM. Adding a fourth ECAD-paired target concentrates the acinar
   autofluorescence bias (see below) on the entire first pass, and leaves no counterpart diversity
   to cross-check against. `panCK` is also acinar-weighted, so it does not escape that bias — but
   it does keep one target off the ECAD web.

---

## 3. Biases carried forward

**Acinar autofluorescence.** Zymogen granules and lipofuscin make exocrine the most
autofluorescent compartment, so an ECAD-high (or panCK-high) negative control **overestimates**
background and every threshold paired to it runs **conservative** — dim positives are lost, false
positives are not gained. The islet-ROI mitigation does not apply: this is a whole-section
compartment pass, so the negative is acinar-dominated by construction. Accept conservative
first-pass cuts and refine at subtyping.

**Prevalence.** All five gates have real prevalence in every section, so the near-zero-prevalence
trap that broke KRT14/TP63/SOX2 does not apply. Neural is sparsest but reliably present.

---

## Configuration as implemented

```python
COMPARTMENT_ORDER = ["Immune", "Epithelial", "Endothelial", "Neural", "Mesenchymal"]  # -> Other

FIRST_PASS_PAIRS = [
    ["ECAD",   "CD31"],
    ["CD31",   "ECAD"],
    ["CD45",   "ECAD"],
    ["PGP9.5", "panCK"],   # substituted for epCAM
    ["VIM",    "ECAD"],    # NOT CD31 — endothelium is vimentin-positive
]
ENDOCRINE_SPLIT_PAIR = ["CHGA", "CD31"]   # inside the Epithelial compartment
```

Two order constraints are enforced as executable assertions in `codex_lineage.validate_order()`,
not left to comment discipline:

- **Immune before Endothelial** — CD31/PECAM1 is on monocytes, neutrophils, NK and naive T cells.
- **Epithelial before Neural** — islet endocrine cells are PGP9.5+/B3TUBB+ and must be cleared first.

`COMPARTMENT_ORDER_EPITHELIAL_FIRST` is provided and validates: ECAD and CD45 swap freely, so
clearing the 85–90% epithelial mass first to shrink every later gate's pool is equally correct.
