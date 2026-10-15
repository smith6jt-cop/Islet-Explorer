# All-donor local delivery

## Phase-2 contract — 2026-09-23

The [current phase-2 contract](donor_expansion_phase2.md) is authoritative.
The six additions are staged separately; original 15-donor IDs, geometry,
measurements, missingness and phenotype labels remain fixed. A validated common
21-donor refit may change derived trajectory, UMAP and Leiden outputs. Earlier
all-21 census, scientific-method and service statements
below describe historical experiments, not the current default app or an
accepted expanded release. Raw staging alone does not make an app-ready dataset.

## Restoration status — 2026-09-23

The local default app is restored at `6b7884b`; the earlier candidate delivery
record below is historical, not current scientific acceptance. Exact original
source-record parity is required before app integration; a future common refit
may change derived trajectory, UMAP and Leiden outputs. The
regenerated schema-4 candidate fails that requirement and needs revision and
evaluation. Its historical restart commands now additionally require
`--baseline-root` and a passing preservation check. No public service cutover
or promotion was performed during the restoration.

## Historical schema-4 candidate — 2026-09-22

At the September 22 checkpoint, the shared app at `http://10.15.152.7:8080/islet-explorer/` served
`data/development/all-donor-20260922-trajectory-fit-support-v1/`. Nginx on port 8080 is the only
network-facing route and proxies to the loopback service on
`127.0.0.1:3840`. The service record contains one-line restart and stop
commands, log paths, the exact app snapshot and the public URL.

The candidate contains seven ND, seven autoantibody-positive and seven T1D
donors, 11,944 retained islet annotations, 3,619,534 canonical raw cell
objects and all 31 shared markers. It preserves 3,592,722 archival source
objects. All 248 wholly unmeasured canonical objects remain present and are
excluded from marker-dependent analyses with recorded reasons; 439 corrected
objects are unavailable for complete-marker analyses. HPAP-078 retains HPAP
provenance and its corrected endocrine-derived boundaries.

Trajectory and Spatial eligibility is recorded per islet and has exactly two
requirements: the endocrine pixel-classifier core does not contact an image
edge, and at least one canonical cell centroid lies inside that classifier
core. There are 11,136 eligible islets across all 21 donors; the smallest donor
contributes 137. Of the 808 ineligible classifier regions, all are empty and
none touch an image edge. Their images and raw measurements remain available.
Downstream phenotype or typing support cannot revoke classifier eligibility.

The saved trajectory is downstream of cell typing and uses core
`log1p(validated base-type count / canonical classifier-core cell count)`
features. The denominator counts every canonical cell centroid in the
classifier-determined endocrine core exactly once.
Disease labels, clinical outcomes, section area, whole-islet marker averages
and earlier pseudotime do not enter the fit. The graph-training subset requires
valid geometry, complete typed-composition features and at least ten measured
core cells; this selects 1,514 islets without changing classifier eligibility.
Twenty deterministic landmarks per donor produce one connected landmark
component. Of 11,136 eligible islets, 7,493 have complete typed-composition
features and 7,490 receive supported primary coordinates after Nyström
projection. The prespecified checks give median leave-one-donor-out absolute
Spearman 0.754 and 67.3% total classifier-eligible coverage, with 21 donor
omissions and 200 bootstrap reconstructions complete.
The ordering is therefore **unstable**. Parameters were not changed to recover
disease ordering; orientation remains arbitrary and is not elapsed disease time.
The raw coordinate contains 1,729 unique values and a 2,234-islet exact tie. The
default 0–1 reference rank preserves raw ordering and spreads exact ties by
SHA-256 of islet ID, giving all 7,490 supported islets a unique fixed midpoint
position. Order within a tie is only for display and descriptive binning;
inferential models retain the raw tied coordinate.

Spatial embeds the 7,493 feature-supported eligible islets from the same typed,
classifier-core-normalized composition matrix. It displays UMAP coordinates, Leiden
resolutions 0.3/0.5/0.8/1.0, donor-status UMAP and cluster composition. The
Trajectory tab displays one point for each of the 7,490 supported primary-axis
islets, nonnegative log1p ratios, the pseudotime strip, restored display
controls and a donor-weighted 31-marker typed-expression heatmap. The default
beta/core view contains 720 nonzero islets representing 6,610 validated beta
cells. Schema-4
drilldown joins donor Parquet cells to exact 20-µm memberships and typed
identities; browser acceptance downloaded 147 cells for a real test islet.

Research results are excluded from the app runtime. Historical result tables
from earlier geometry are not current-candidate evidence. The independent
schema-4 validator and browser regression pass, but `scientific_acceptance`
and UFRC deployment acceptance remain false.

## Historical delivery records

**Previous revision: 2026-09-17.** The earlier preview on 3839 used
`data/development/all-donor-20260917/` and schema 3. Research results are offline,
all 31 heatmap markers are enabled, and the newly identified islet-identity
problem is visible in the app. Ratios still require validated parent gates.
See the [current implementation, review evidence and restart record](trajectory_revision_20260917.md).
The dated 2026-09-16 delivery and earlier repair record below are preserved as
history; their Research tab, schema, capacity and passing evidence describe that
earlier candidate and do not supersede the current record.

The separate candidate is `data/development/all-donor-20260916/`. Its local
preview is **http://127.0.0.1:3839/islet-explorer/**. Port 3838 was occupied;
its service was not stopped or reconfigured. The candidate uses its own Shiny
Server, application snapshot, static TIFF routes and logs. This is local
delivery, with scientific acceptance and UFRC deployment acceptance both false.
The frozen 15-donor `5809be4` release and its dated evidence remain intact.

## Cohort and measurements

The [explicit cohort manifest](../scripts/research/cohort-v1.csv) selects 21
donors: seven ND, seven autoantibody-positive and seven T1D. HPAP-078 retains
HPAP provenance (analysis ID `78`, viewer ID `0078`); the other 20 donors are
nPOD sources. Corrected 6539 comes from the recovered acquisition, excluding
the rejected duplicate. Alternate CX_PAN-6558 and 6520-01B remain excluded.

| Candidate quantity | Count |
|---|---:|
| Raw cell objects retained | 3,592,722 |
| Wholly unmeasured objects retained as NaN | 201 |
| Shared raw markers | 31 |
| Islet regions and available cell drill-down files | 8,371 |
| Islets eligible for the primary trajectory | 3,499 |
| Compartment/band measurement rows | 50,226 |

The 201 unmeasured objects are excluded from marker-dependent analyses with
recorded reasons. Raw intensities, archival labels and images are retained.
Measurements use 0.3774 µm/pixel, tissue clipping and 10/20/30-µm peri-islet
bands that exclude neighboring cores. Structural density is positive mask
area divided by valid compartment area. Primary spatial tests exclude
overlapping 30-µm footprints and 6539's still-unreviewed exclusion annotations.

The current retirement policy applies to nine obsolete labels. Validation
audits each retained label's complete parent-gate chain; reconstructing all
21 historical labels is no longer a prerequisite. The available upstream
gates do not certify subtype assignments. Consequently the alpha/beta/delta
union denominator and dependent trajectory ratios remain explicitly unavailable.
Zero denominators remain missing. The app retains immune and structural
cell-type choices, core/peri/both numerator choices, raw marker heatmaps and
cell-level marker downloads. No retired Endocrine denominator is substituted.

## Fixed trajectory

The offline method uses core and 20-µm peri-islet median `log1p` raw intensities,
excluding INS/GCG/SST/CHGA. Clinical outcomes, endocrine ratios, disease labels
and prior coordinates do not enter the fit or endpoint choice. Landmark
sampling is equal per donor across section-area quartiles, capped at 200 and
limited by the smallest eligible donor set. Scaling, cosine k=15 neighbors,
ten diffusion components, seed 20260916, two-sweep endpoints and hash orientation
are fixed. Non-landmarks use Nyström projection with recorded support.
The smallest eligible donor contributes eight islets, so the fitted landmark
set contains eight per donor, 168 in total.

All 3,499 eligible islets receive supported primary-component coordinates;
the fitted graph has one component. Median leave-one-donor-out absolute
Spearman correlation is **0.9143**, passing the prespecified 0.8 threshold.
All 21 donor exclusions, 200 donor bootstraps and five parameter/seed
sensitivities completed. Bootstrap median absolute rank correlation is
approximately **0.759**, with a wide range; passing the specified thresholds
does not establish biological robustness or absence of bias. Donor associations
are strong, and one acquisition per donor prevents separating acquisition
from donor effects. Orientation is arbitrary and is not elapsed disease time.

The app reads saved features, transformations, landmarks, graph, root and
coordinates. Donor selection never refits or reranks the axis. Its downloads,
heatmap and highlighted points share the selected donor; the full graph remains
faintly visible. Descriptive curves have no population confidence intervals.
No inferential pseudotime outcome associations have been published. Such tests
still require outcome-module/gating exclusion, preprocessing and projection
within donor folds, trajectory-reconstructing bootstrap uncertainty and their
own multiple-testing family.

## Primary results and follow-ups

The [versioned registry](../scripts/research/hypotheses-v1.json) specifies ten
primary tests, T1D versus ND, donor replication and BH adjustment across all
ten slots. The Research results tab provides effect estimates, bootstrap
intervals, p/q values, donor/islet counts, measurement downloads, plots and
representative-islet links. Estimates remain pending scientific review.

### Why the primary models include 12–13 donors

The all-donor dataset and fixed trajectory include all **21 donors**. The
primary hypothesis contrast is **T1D versus ND**, starting with seven T1D and
seven ND donors (14); the seven Aab-positive donors are reserved for secondary
contrasts, which have not yet been run.

ND donor **6539** is excluded from these primary endpoint tables because
the corrected acquisition's exclusion/boundary annotation review is pending.
Its image, raw measurements and trajectory coordinates are retained. The
remaining five estimated primary models (H1/H2/H3/H5/H10) therefore use
**six ND plus seven T1D donors = 13**.

H9 additionally excludes T1D donor **6520**. Its implementation requires at
least five measured INS/GCG/SST-positive cells per islet, exact agreement with
the exported hormone-positive union count, and at least three qualifying
islets per donor. All 51 assessable 6520 cores have only 0–4 such cells;
all 51 union counts agree with the source export. Thus none qualifies for
H9, leaving **six ND plus six T1D donors = 12**. Donor 6520 remains in the
other five estimated primary tests. H9's comparison is conditional on having
enough measurable endocrine cells; this exclusion must remain visible when
interpreting heterogeneity.

The donor-level audit is saved separately from the immutable result artifacts
in `evidence/hypothesis_donor_inclusion.csv` (verified 2026-09-17).

| Test | Donors / islets | Effect | BH q | Status |
|---|---:|---:|---:|---|
| H1 size × insulin composition | 13 / 810 | −0.5304 log-odds interaction | 0.3846 | Estimated |
| H2 nerve compartment contrast | 13 / 955 | −0.02288 area-fraction interaction | 0.06517 | Estimated |
| H3 capillary compartment contrast | 13 / 955 | −0.009876 area-fraction interaction | 0.06517 | Estimated |
| H4 lymphatic–immune association | — | — | — | Validated CD45/PDPN gates unavailable |
| H5 neurovascular coupling | 13 / 890 | 0.04849 excess-proximity contrast | 0.3846 | Estimated; uniform-location spatial null |
| H6 cytotoxic-cell access | — | — | — | Validated CD3e/CD8a parent gates unavailable |
| H7 basement membrane | — | — | — | Validated ColIV boundary mask/CD45 gate unavailable |
| H8 ductal niches | — | — | — | Validated CK19/panCK duct mask unavailable |
| H9 endocrine heterogeneity | 12 / 227 | 0.3509 log-dispersion contrast | 0.5173 | Estimated using all 31 markers |
| H10 added predictive information | 13 / 425 | −0.2262 balanced-accuracy change | 1.000 | Estimated with nested donor-held-out fitting |

None passes BH q<0.05. H1's operational hormone-positive union comes from the
QuPath INS/GCG/SST export and does not certify phenotype subtypes. H9 additionally
requires exact reconstructed union-count agreement and adequate measured
endocrine-cell counts. H10 repeats its entire nested fitting procedure for
200 donor-label permutations and 200 grouped donor bootstrap attempts; baseline
balanced accuracy is 1.000 and full-profile accuracy is 0.7738 in these 13 donors.
Small samples and donor influence limit interpretation.

F1 follows the prespecified unsupported-primary branch: H1 adjusted for both
structural densities (effect −0.5504, bootstrap CI −1.2056 to −0.0561,
p=0.09059; exploratory). F2/F3 await missing primary results; unavailable tests
are not interpreted as evidence of no association. F4 selects the prespecified
size-plus-structural-density fallback, but testing against INS-positive tissue
burden awaits a validated tissue-area measurement. Cell counts do not replace it.

Scientific work still needed includes retained-gate and mask review, corrected
6539 annotation review, shape-preserving spatial-null and remaining band/
segmentation sensitivities, and the unavailable endpoints/follow-ups. Aab-positive
and clinical comparisons are secondary and currently untested; preserve missing
and censored clinical values. Full implementation details and reproduction are
in the [analysis README](../scripts/research/README.md).

## Runtime and provenance

`manifest.json` declares runtime schema **2**, full artifact SHA-256 hashes,
source hashes, 21 images and explicit unavailable analyses. It adds cohort,
trajectory and research artifacts; raw cell drill-downs are required. The R
reader verifies startup artifacts independently, using portable byte ordering
for schema 2. Schema 1 retains its original digest and artifact conventions.
Deployment verification and the production R profile reject this candidate
while `scientific_acceptance` is false.
The final manifest SHA-256 is
`4826e871cd3ac6517ed6a78ea1da8bcf3f47949f77e4fedaae61def9b96a7bf5`.

| Measured storage | Bytes | GiB |
|---|---:|---:|
| Candidate images and sidecars | 163,098,139,438 | 151.897 |
| Complete declared runtime | 164,130,834,281 | 152.859 |
| Local candidate, including analysis/source material | 170,963,414,830 | 159.222 |
| Frozen original plus candidate runtime | 281,223,614,814 | 261.910 |
| Frozen original plus candidate and one staging copy | 445,354,449,095 | 414.769 |

Local logical sizes count hardlinked files; portable copies require the full
declared bytes. Never modify hardlinked raw inputs or images. Geometry occupies
786,821,120 bytes and per-islet cell CSVs 216,569,294 bytes. The two-runtime-copy
capacity field excludes the separately retained original; use the last row
above when planning migration. Allow additional OCI, logs, temporary/spill and
growth space. The requested 1-TiB allocation has not been confirmed or load-tested.

The compartment-building stages recorded 259.6 seconds cumulatively on this
workstation, excluding upstream image processing, hypothesis fitting, DPT
validation and hashing. Full image hashing was disk-bound and took more than
ten minutes. The final manifest-backed worker's observed RSS was 765,920 KiB
(0.730 GiB); earlier preview observations reached approximately 1 GiB. The final
browser connected in 3.231 seconds and completed its checks in 17.211 seconds.
These are workstation observations, not peak-memory measurements, a five-session
load test or UFRC runtime acceptance.

## Local operation and verification

From the repository root, restart only this candidate with:

```bash
python3 scripts/research/serve_preview.py \
  --data-root data/development/all-donor-20260916 \
  --app-root app/shiny_app \
  --python /home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  --restart
```

Read `preview/service.json` for the actual URL, current PID, restart/stop argv,
application snapshot and log paths. Logs are `preview/server.log` and
`preview/logs/`. If 3838 is occupied, the launcher uses the first free loopback
port from 3839. Stopping this candidate leaves the existing service alone.
The snapshot can be rebuilt from this checkout; raw/derived candidate artifacts
remain separate. Production rollback continues to use the frozen release and
the [deployment rollback procedure](../deploy/README.md).

Independent candidate validation checks all donors/raw marker missingness,
cohort joins, structural units and fraction bounds, all eligible cell drill-down
paths, exclusion policy, feature coverage and primary result availability.
Python regression tests cover deterministic fits, label/ratio invariance,
disconnected components, projection alignment, outcome-gate exclusion, grouped
donor folds and spatial overlap handling. R tests cover fixed donor coordinates,
zero-denominator ratios and both runtime schemas. Browser evidence checks the
21 donors, donor-filtered downloads, immutable coordinates, raw-cell drill-down,
HPAP image click-through, corrected-6539 HTTP Range, results and accessibility.

The completed local checks passed **47 Python tests**, **204 R assertions in
45 tests**, and synthetic 21-image activation/15-image rollback. Browser
checks use Node 22.14.0 and Chrome 153.0.8010.47, with no axe violations across
the five app tabs. These checks cover the candidate's explicit unavailable
states and do not certify the missing biological inputs.

The restricted candidate `evidence/` directory contains the final checks,
manifest identity, source inventory and browser versions. No frozen passing
record is reused as candidate validation. Automated accessibility checks do
not replace manual screen-reader, zoom or scientific-plot review. A clean
release revision, matching OCI image/SBOM/scan, production source parity,
UFRC host/auth/proxy verification, concurrency measurements and institutional
acceptance remain outstanding. Existing UFRC request/assessment identifiers
and ownership requirements remain in the [handoff](../PUBAPPS_UFRC_HANDOFF.md).

Key machine-readable records are `evidence/local-acceptance.json`,
`evidence/source-and-manifest-identity.json`, `candidate_validation.json`,
`capacity.json`, `trajectory/validation.json` and `research/results.csv`.
The checkout remains uncommitted; the evidence records the base revision and
exact runtime source hashes rather than claiming a clean release revision.

## Trajectory UI repair — 2026-09-17

The expected unavailable-ratio condition was displayed through Shiny's plot
error mechanism. It now appears as a compact availability notice, keeping the
fixed diffusion graph and marker heatmap visible. This presentation repair
does not certify phenotype gates or create the unavailable ratio results.
A separate zero-row assignment error when filters selected no islets was also
fixed. Empty downloads retain their headers, and empty heatmap/marker selections
show informational prompts and recover when the selection is restored.

Fresh verification passed 225 R assertions across 48 tests, runtime lint,
and browser checks including empty-filter recovery, fixed coordinates,
downloads, image links and zero axe violations on the five tabs. The app
was restarted at the same local URL. Scientific artifacts and their manifest
were unchanged. Dated evidence and runtime source hashes are retained in
`evidence/trajectory-fix-20260917/`; the original delivery evidence remains intact.
## Candidate delivery record — 2026-09-22

The schema-4 all-donor candidate is the only network-facing app at
`http://10.15.152.7:8080/islet-explorer/`. Nginx proxies that URL to the
loopback backend on `127.0.0.1:3840`; the older 3838 and 3839 services are
stopped and disabled. The service record, logs, and nginx configuration backup
are recorded in `data/development/all-donor-20260922-trajectory-fit-support-v1/preview/service.json`,
`data/development/all-donor-20260922-trajectory-fit-support-v1/preview/server.log`, and
`/etc/nginx/conf.d/islet-explorer.conf.before-single-20260921`.

Restart the backend by copying this one line into a shell:

```bash
cd /data/panc_CODEX/Islet-Explorer && /home/smith6jt/miniconda3/envs/scvi-env/bin/python3.13 scripts/research/serve_preview.py --data-root data/development/all-donor-20260922-trajectory-fit-support-v1 --app-root app/shiny_app --python /home/smith6jt/miniconda3/envs/scvi-env/bin/python --port 3840 --restart
```

The public URL remains unchanged after a backend restart. If nginx itself is
restarted, reload it with `pkexec systemctl reload nginx`.

The preview and public proxy respond with HTTP 200. Browser acceptance verifies
the Spatial UMAP and Leiden views, the 31-marker heatmap, fixed donor filtering,
ratio arithmetic, islet click-through, cell downloads and the absence of a
Research tab. The current fixed trajectory records all 21 donors, 11,136
classifier-eligible islets, 7,493 feature-supported rows and 7,490 supported
primary coordinates. Research results remain offline and are not part of the
runtime. Machine-readable browser
evidence and its screenshot are in
`evidence/browser-schema4-trajectory.{json,png}` under the candidate root.
