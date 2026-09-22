# Trajectory and islet review revision — 2026-09-17

The revised local app is **http://127.0.0.1:3839/islet-explorer/**, backed by
`data/development/all-donor-20260917/`. All 21 donors remain discoverable. The
existing service on 3838 and frozen 15-donor release remain intact. Research
results have been removed from the app, and all 31 heatmap markers are selected
by default. **The trajectory is not scientifically accepted.** Neither the
islet-identity problem nor the missing validated endocrine denominator is solved
by the working heatmap and graph.

## Islet-identity finding

The owner's QuPath review identified many similarly sized HPAP-078 annotations
without INS/GCG/SST expression. The audit found **348 of 1,234 HPAP candidates**
with no cells above the literal exported hormone thresholds. **143** were in
the previous trajectory-eligible set. Across all 21 donors, **1,755 of 8,371**
candidates have this flag, including **508** previously eligible regions.
These are review flags, not a determination that all hormone-negative regions
are non-islets.

Six examples of the repeated small outlines are Islet_130, Islet_206,
Islet_1148, Islet_124, Islet_161 and Islet_179. Their bounding boxes are roughly
25–32 µm wide. Fixed-display raw-channel crops show little hormone signal in
these examples; CHGA and tissue context still need inspection. The audit follows
each saved mask to its classifier candidate and SAM prompt. All **2,989 of
2,989** masks from the six additions match their SAM exports, including all
1,234 HPAP masks. The failure is not an image-coordinate offset during import.

The classifier proposes candidate regions; their padded bounding boxes prompt
SAM on a downsampled INS/GCG/SST rendering. SAM's six-pixel grid partly explains
repeated dimensions. The historical importer checked minimal mask geometry and
assigned the class `Islet`, without confirming biological identity. A completed
segmentation run and its mask-confidence score did not validate those regions.

The existing fit and numerical stability checks are preserved for inspection.
`trajectory/quality_review.json` now takes precedence in the app's status:
islet identity requires review. Future command-line fits require per-islet
review decisions bound to the source geometry hashes. Confirmed identity
corrections require a new offline fit and validation. Geometry edits additionally
require rebuilding cell membership, bands, marker summaries and dependent
research analyses. Display filters never initiate a refit.

Review evidence is outside the app:

- `data/reviews/islet-identity-20260917/`: all-donor flags, donor summaries,
  `islet_review_decisions.csv`, and donor-specific flagged GeoJSON overlays.
- `.../sam-shape-review/78/small_candidate_raw_channels.pdf` and `.png`:
  the six small examples, with saved outlines on raw DAPI/INS/GCG/SST/CHGA.
- `.../sam-shape-review/sam_geometry_alignment_corrected.csv`: source-mask
  alignment for all six additions. This corrects the first audit's strict
  vertex comparison, which failed to account for QuPath removing collinear
  vertices. The first report remains preserved; it is not evidence of misalignment.
- `data/reviews/analysis-regions-20260917/`: exact tissue-clipped cores and
  20-µm peri-islet bands, excluding neighboring cores, for all 21 donors.
  Every exported compartment area is checked against the count artifact.

The editable native project remains
`data/reviews/qupath-six-donors-20260917/project.qpproj`. This revision did not
modify its object files or the owner's edits. Import optional review overlays
only into the matching image; do not import duplicates repeatedly. Record
`confirmed_islet`, `non_islet`, `uncertain` or `unreviewed` decisions and the
review basis. Absence of hormone signal alone is not an automatic deletion rule.

## Count engine and remaining gate decisions

The evidence registry now covers all **21 × 31** donor/marker combinations.
Application of reviewed retained rules no longer requires an archival phenotype
label. Archival-label reconciliation is retained as a separate audit. Complete
parent dependencies, finite thresholds and biological validation evidence are
required; unsupported assignments and counts remain missing.

No donor currently has a certified alpha/beta/delta union denominator under
the retained rules: **0 of 252** donor/retained-type pairs are validated. The
endocrine parent chain needs INS/GCG/SST plus CD20, CD3e, CD68, CD163, CD8a, CD4
and HLADR evidence. `data/reviews/phenotype-gates-20260917/` contains the registry,
per-donor intensity plots, nullable cell assignments, archival audits and review
templates. `endocrine_gate_review_decisions.csv` prioritizes the **165** pending
donor/marker decisions needed for the endocrine parent chain. Plots and exported
thresholds support review but do not constitute biological validation.

The app reads dedicated numerator and endocrine-denominator artifacts, exports
their validation status and availability reasons, and retains all 12 permitted
cell types. Core/peri/both numerators share the unique core endocrine union.
Validated zero numerators stay zero; zero or unknown denominators stay missing.
The 201 wholly unmeasured raw objects remain archived and excluded from marker
analyses. All 8,371 core cell memberships agree with the source membership cache.

## App and offline reports

The trajectory heatmap offers all 31 markers, including INS/GCG/SST/CHGA, with
Select all and Clear controls. Rows preserve panel order; the plotted measure
is the mean of per-islet median log1p core intensity within fixed trajectory
bins. Donor filters apply to heatmaps, downloads and graph highlights while
coordinates, reference percentiles and axis limits remain fixed.

Research reports are at
`data/reports/all-donor-research/20260917/index.html`. H1–H10 folders contain
available estimates, intervals, p/q values, counts, charts in PNG/PDF and source
image references; unavailable tests have explicit reasons. These are the
previously computed estimates, copied without refitting and now explicitly
pending islet-identity review. They must not be interpreted as accepted results.
All offline research runners require an explicit report root outside app data.

Runtime schema **3** requires count/availability and marker-panel artifacts,
and removes the research-data dependency. Schema 1/2 readers remain compatible
with preserved releases. There is no Research tab, result download, result
resource route or research folder in this candidate. Both the production R
profile and deployment verifier continue to reject scientific acceptance=false.

## Fresh verification and capacity

The revision passed 50 R tests with 241 assertions, 54 Python tests (20 research
contracts plus 34 cohort/metadata tests), runtime lint, both synthetic deployment
test suites, candidate validation, independent R startup hash verification and browser
checks. The browser checked 21 donors/images, all 31 marker defaults/exports,
donor-filtered downloads, fixed coordinates, empty-state recovery, raw-cell
download, HPAP image click-through, corrected-6539 TIFF Range, absence of research
UI/routes and zero axe violations across the four tabs. **Actual populated
ratio curves remain unverified because their scientific inputs are unavailable.**
Numerical ratio behavior is covered by synthetic tests; an all-missing notice
does not count as successful scientific unblocking.

Nine fitted trajectory artifacts are byte-identical to the preceding candidate.
The new full runtime manifest contains 13 artifacts and has SHA-256
`f18477f37dc6bafc9e4f044dc867b49fb1329dab15caadd8b985598d93f9aedf`.
Fresh evidence is under `data/development/all-donor-20260917/evidence/`, with
logs in `logs/` and `preview/logs/`. Earlier passing evidence remains dated
to the earlier candidate and is not reused as this revision's acceptance.

| Measured quantity | Bytes / observation |
|---|---:|
| Declared runtime | 164,147,382,976 bytes; 152.874 GiB |
| Images and sidecars | 163,098,139,438 bytes; 151.897 GiB |
| Candidate logical size at finalization | 170,975,613,581 bytes; 159.233 GiB |
| Frozen original plus candidate | 281,240,163,509 bytes; 261.925 GiB |
| Frozen original plus candidate plus one staging copy | 445,387,546,485 bytes; 414.799 GiB |
| Observed new preview R worker RSS after browser checks | 736,172 KiB; about 0.70 GiB |
| Browser connection / entire run | 4.595 seconds / 27.777 seconds |

Logical sizes count immutable hardlinks. Revised derived files are independent
copies. Review projects/reports, application snapshots, growing logs, OCI images
and spill/growth space require additional capacity. These observations are not
peak-memory, concurrency, assigned-host or UFRC acceptance measurements. The
259.6-second measurement-builder timing retained in `capacity.json` belongs to
the original measurement build; it is not a new full-pipeline benchmark.

## Local restart and rollback

```bash
python3 scripts/research/serve_preview.py \
  --data-root data/development/all-donor-20260917 \
  --app-root app/shiny_app \
  --python /home/smith6jt/miniconda3/envs/scvi-env/bin/python --restart
```

`preview/service.json` records the actual URL, PID, restart/stop argv and logs.
The launcher checks the exact recorded configuration before stopping a preview.
It replaced only the earlier candidate on 3839; the service on 3838 was preserved.

For an exact local rollback, stop only the current recorded preview and start
`/opt/shiny-server/bin/shiny-server` with the preserved configuration at
`data/development/all-donor-20260916/preview/shiny-server.conf`; its application
snapshot includes the original isolated environment profile. Check the recorded
PID/configuration before stopping a process and verify the old port is free.
Do not use the old data with a rebuilt snapshot when claiming exact code rollback.
The frozen production release and its deployment rollback instructions are
unchanged. Local delivery does not establish scientific or UFRC acceptance.
