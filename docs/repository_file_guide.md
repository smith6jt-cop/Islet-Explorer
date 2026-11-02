# Islet Explorer: file-by-file reference

**Current additions, 2026-11-02:** [the workflow user guide](workflow_user_guide.md) covers
the present source-to-app chain. `scripts/uniform/typing/calibration_context.py` pins reference
context; `typing/structure_update.py` pins a structure revision's calibration and baseline;
`scripts/uniform/audit_structure_update.py` verifies every component after identity generation.
The [component comparison](v19_v20_component_comparison_20261102.md) and
[cross-repository change notice](typing_change_notice_20261102.md) describe the selected sources,
fixes, evidence and remaining work. The file count below belongs to its stated September revision.

Inspected 2026-09-10 at Git `000e970cac91220f71200b5240f4729a3d32d44e`. This catalog explicitly describes all **225 tracked entries** at that revision, including the submodule reference. It supplements [the maintainer walkthrough](maintainer_walkthrough.md).

The subsequent offline expansion workflow and its separate data staging are
documented in [the expansion status](cohort_expansion_status.md) and
[processing README](../scripts/cohort_expansion/README.md). This catalog's file
count and frozen-data examples retain their stated revision.

The local searchable TSV inventory under `tmp/maintainer-training-20260910/` also includes ignored data, installed dependencies, temporary releases, and symlinks. It records metadata only. Git object-store internals and managed assistant directories are represented as boundaries rather than enumerated internally. Directory symlinks are recorded without following them.

Runtime inclusion is controlled by the Containerfile/build-context whitelist and the release-data builder, not simply by Git tracking. Some tests/docs are packaged under the app tree without being executed. A lifecycle label identifies purpose rather than certifying exact image contents.

A separate review register records candidate assistant references and AI-feature/control references with line numbers. It is a lexical aid, not a decision to edit. Your request makes removal conditional on your file review, and the precise meaning of “AI attributes” is awaiting clarification. Microscopy **CODEX** and the browser `navigator` API are not AI authorship markers.

**Repository root and tooling**

| File | Role |
|---|---|
| [.claude/hooks/session-start.sh](../.claude/hooks/session-start.sh) | Coding-assistant session initialization hook; development tooling. |
| [.claude/settings.json](../.claude/settings.json) | Shared coding-assistant hook/tool configuration. |
| [.github/copilot-instructions.md](../.github/copilot-instructions.md) | Repository guidance for a coding assistant. |
| [.github/workflows/pubapps-release.yml](../.github/workflows/pubapps-release.yml) | Builds and scans the immutable OCI image, exports it, and retains evidence; does not deploy the VM. |
| [.github/workflows/shiny-app-ci.yml](../.github/workflows/shiny-app-ci.yml) | Source parsing, R tests/lint, Python analysis, browser SAST, and deployment contract CI. |
| [.gitignore](../.gitignore) | Excludes bulk data, local settings, environments, caches, and staging from Git; explicitly includes selected provenance/metadata files. |
| [.gitmodules](../.gitmodules) | Defines the separately versioned Phenocycler_Analysis submodule and its upstream URL/branch. |
| [CLAUDE.md](../CLAUDE.md) | Broad project development instructions, scientific invariants, history, and links; contains assistant-specific framing. |
| [Islet-Explorer.code-workspace](../Islet-Explorer.code-workspace) | VS Code workspace configuration. |
| [LICENSE](../LICENSE) | Software license; does not authorize redistribution of research data. |
| [PUBAPPS_UFRC_HANDOFF.md](../PUBAPPS_UFRC_HANDOFF.md) | Request facts, intended access/data profile, owner responsibilities, copy-ready form material, and go-live checklist. |
| [README.md](../README.md) | Project overview, features, development entrypoint, data architecture, and links to current production instructions. |
| [UFIT_PubApps_Instance_Request.docx](../UFIT_PubApps_Instance_Request.docx) | Populated UFRC request form plus supplied policy, project plan, and development standards; acknowledgment still needs the responsible person's review. |

**Application entrypoint and shared runtime**

| File | Role |
|---|---|
| [app/docs/logo.png](../app/docs/logo.png) | Application logo source asset; served through the www/logo.png link. |
| [app/shiny_app/CLAUDE.md](../app/shiny_app/CLAUDE.md) | Detailed application development notes and scientific/UI invariants with assistant-specific framing. |
| [app/shiny_app/Channel_names](../app/shiny_app/Channel_names) | Marker/channel naming reference. |
| [app/shiny_app/R/00_audit.R](../app/shiny_app/R/00_audit.R) | Allowlisted one-line JSON event emission and shared validation primitives; logs fixed codes rather than raw scientific inputs. |
| [app/shiny_app/R/00_globals.R](../app/shiny_app/R/00_globals.R) | Package loading, production debug rejection, data/Python paths, palettes, scientific constants, feature flags, and common UI helpers. |
| [app/shiny_app/R/01_runtime_data_backend.R](../app/shiny_app/R/01_runtime_data_backend.R) | Process-wide cache, embedded DuckDB queries, release schema/path/hash validation, and shared H5AD loading. |
| [app/shiny_app/R/ai_helpers.R](../app/shiny_app/R/ai_helpers.R) | Retained optional AI credential/provider helpers; production disables AI and blocks credential discovery/request paths. |
| [app/shiny_app/R/barlow_pseudotime.R](../app/shiny_app/R/barlow_pseudotime.R) | Composition-based diffusion pseudotime support and caching for the selected population. |
| [app/shiny_app/R/data_loading.R](../app/shiny_app/R/data_loading.R) | Development H5AD/Excel loading, scientific preparation, and production prepared-release loading/cache. |
| [app/shiny_app/R/drilldown_helpers.R](../app/shiny_app/R/drilldown_helpers.R) | Per-islet cell CSV lookup, caching, and cell-level plotting/summaries. |
| [app/shiny_app/R/mod_ai_assistant_server.R](../app/shiny_app/R/mod_ai_assistant_server.R) | AI Assistant (disabled): reactive logic for retained provider-backed assistant interface and request logic. |
| [app/shiny_app/R/mod_ai_assistant_ui.R](../app/shiny_app/R/mod_ai_assistant_ui.R) | AI Assistant (disabled): UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_inflammation_server.R](../app/shiny_app/R/mod_inflammation_server.R) | Inflammation (disabled): reactive logic for retained cohort-level abundance/infiltration analyses. |
| [app/shiny_app/R/mod_inflammation_ui.R](../app/shiny_app/R/mod_inflammation_ui.R) | Inflammation (disabled): UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_plot_server.R](../app/shiny_app/R/mod_plot_server.R) | Plot: reactive logic for filtering, size/marker/composition plots, exports, and selected-islet events. |
| [app/shiny_app/R/mod_plot_ui.R](../app/shiny_app/R/mod_plot_ui.R) | Plot: UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_spatial_server.R](../app/shiny_app/R/mod_spatial_server.R) | Spatial: reactive logic for neighborhood, donor tissue, cell distance/density/composition views, and exports. |
| [app/shiny_app/R/mod_spatial_ui.R](../app/shiny_app/R/mod_spatial_ui.R) | Spatial: UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_statistics_server.R](../app/shiny_app/R/mod_statistics_server.R) | Statistics (disabled): reactive logic for retained statistical comparisons and exports. |
| [app/shiny_app/R/mod_statistics_ui.R](../app/shiny_app/R/mod_statistics_ui.R) | Statistics (disabled): UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_trajectory_server.R](../app/shiny_app/R/mod_trajectory_server.R) | Trajectory: reactive logic for composition trajectory/PAGA, count-range filtering, plots, linked selections, and exports. |
| [app/shiny_app/R/mod_trajectory_ui.R](../app/shiny_app/R/mod_trajectory_ui.R) | Trajectory: UI controls, labels, and output locations. |
| [app/shiny_app/R/mod_viewer_server.R](../app/shiny_app/R/mod_viewer_server.R) | Viewer: reactive logic for approved image selection and local Avivator embedding. |
| [app/shiny_app/R/mod_viewer_ui.R](../app/shiny_app/R/mod_viewer_ui.R) | Viewer: UI controls, labels, and output locations. |
| [app/shiny_app/R/paga_layout.R](../app/shiny_app/R/paga_layout.R) | R reticulate bridge to the Python PAGA/layout module and cached results. |
| [app/shiny_app/R/segmentation_helpers.R](../app/shiny_app/R/segmentation_helpers.R) | Segmentation lookup, production indexed GeoPackage access, geometry/coordinate handling, and plotting support. |
| [app/shiny_app/R/spatial_helpers.R](../app/shiny_app/R/spatial_helpers.R) | Spatial distance, donor-tissue, and related queries; production Parquet/DuckDB paths and development alternatives. |
| [app/shiny_app/R/utils_safe_join.R](../app/shiny_app/R/utils_safe_join.R) | Join helpers preserving the intended scientific identifiers and compatibility behavior. |
| [app/shiny_app/R/utils_stats.R](../app/shiny_app/R/utils_stats.R) | Reusable statistical/summary helpers; individual helper use spans modules, including retained optional analyses. |
| [app/shiny_app/R/viewer_helpers.R](../app/shiny_app/R/viewer_helpers.R) | Local image discovery, safe file/URL handling, and viewer configuration. |
| [app/shiny_app/README.md](../app/shiny_app/README.md) | Older app introduction with stale launch/data-path examples; use root README and deploy/README for current execution. |
| [app/shiny_app/app.R](../app/shiny_app/app.R) | Live Shiny entrypoint: top-level UI, process-wide prepared data, per-session reactives, module wiring, and common drilldown/export handlers. |
| [app/shiny_app/docs/spatial.md](../app/shiny_app/docs/spatial.md) | Scientific and UI explanation for Spatial analyses. |
| [app/shiny_app/docs/statistics.md](../app/shiny_app/docs/statistics.md) | Scientific explanation for the retained, disabled Statistics module. |
| [app/shiny_app/docs/trajectory.md](../app/shiny_app/docs/trajectory.md) | Current trajectory/PAGA methodology, feature construction, count-range semantics, and interpretation. |
| [app/shiny_app/python/paga_pseudotime.py](../app/shiny_app/python/paga_pseudotime.py) | Runtime Python computation of composition neighbors, Leiden/PAGA, diffusion pseudotime, and force-directed layout, called from R. |
| [app/shiny_app/restart.txt](../app/shiny_app/restart.txt) | Legacy Shiny application restart marker; excluded from production image context and not the Quadlet lifecycle control. |

**Application tests**

| File | Role |
|---|---|
| [app/shiny_app/tests/README.md](../app/shiny_app/tests/README.md) | How to run/interpret the app tests; historical assistant/development references are present. |
| [app/shiny_app/tests/manual/stats_audit.R](../app/shiny_app/tests/manual/stats_audit.R) | Manual statistical audit harness for retained analysis behavior. |
| [app/shiny_app/tests/parse_check.R](../app/shiny_app/tests/parse_check.R) | Parses R application sources without launching a reviewer session. |
| [app/shiny_app/tests/run_tests.R](../app/shiny_app/tests/run_tests.R) | Application testthat suite runner. |
| [app/shiny_app/tests/testthat/helper-setup.R](../app/shiny_app/tests/testthat/helper-setup.R) | Shared test fixtures and helper loading. |
| [app/shiny_app/tests/testthat/test-application-logging.R](../app/shiny_app/tests/testthat/test-application-logging.R) | Checks fixed-field JSON logs and rejection of unsafe log inputs. |
| [app/shiny_app/tests/testthat/test-runtime-data-backend.R](../app/shiny_app/tests/testthat/test-runtime-data-backend.R) | Checks production manifest/data backend behavior and development/production separation. |
| [app/shiny_app/tests/testthat/test-security-boundaries.R](../app/shiny_app/tests/testthat/test-security-boundaries.R) | Checks selected input/path/production security boundaries. |
| [app/shiny_app/tests/testthat/test-trajectory-helpers.R](../app/shiny_app/tests/testthat/test-trajectory-helpers.R) | Checks trajectory helper computation/contracts. |
| [app/shiny_app/tests/testthat/test-utils_stats.R](../app/shiny_app/tests/testthat/test-utils_stats.R) | Checks statistical helpers against expected outcomes. |
| [app/shiny_app/tests/testthat/test-viewer-readiness.R](../app/shiny_app/tests/testthat/test-viewer-readiness.R) | Checks viewer readiness and controlled behavior when image inputs are missing or invalid. |

**Browser assets**

| File | Role |
|---|---|
| [app/shiny_app/www/apple-touch-icon.png](../app/shiny_app/www/apple-touch-icon.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/avivator/VENDOR_PROVENANCE.md](../app/shiny_app/www/avivator/VENDOR_PROVENANCE.md) | Provenance evidence and explicit limits of upstream revision/dependency reconstruction for the checked-in vendor bundle. |
| [app/shiny_app/www/avivator/VENDOR_SHA256SUMS](../app/shiny_app/www/avivator/VENDOR_SHA256SUMS) | Expected SHA-256 hashes of retained vendor assets. |
| [app/shiny_app/www/avivator/assets/deflate-Bj5r1ZIu.js](../app/shiny_app/www/avivator/assets/deflate-Bj5r1ZIu.js) | Vendored Deflate image decoder chunk. |
| [app/shiny_app/www/avivator/assets/index-DBmmXfNW.js](../app/shiny_app/www/avivator/assets/index-DBmmXfNW.js) | Prebuilt third-party viewer application bundle; direct modifications affect pinned vendor hash/provenance. |
| [app/shiny_app/www/avivator/assets/index-jOqUmYCe.css](../app/shiny_app/www/avivator/assets/index-jOqUmYCe.css) | Prebuilt viewer stylesheet. |
| [app/shiny_app/www/avivator/assets/jpeg-BXvFeBhq.js](../app/shiny_app/www/avivator/assets/jpeg-BXvFeBhq.js) | Vendored JPEG image decoder chunk. |
| [app/shiny_app/www/avivator/assets/lerc-DeOxovhr.js](../app/shiny_app/www/avivator/assets/lerc-DeOxovhr.js) | Vendored LERC image decoder chunk. |
| [app/shiny_app/www/avivator/assets/lzw-4WQwzCSE.js](../app/shiny_app/www/avivator/assets/lzw-4WQwzCSE.js) | Vendored LZW image decoder chunk. |
| [app/shiny_app/www/avivator/assets/packbits-D01fq4dm.js](../app/shiny_app/www/avivator/assets/packbits-D01fq4dm.js) | Vendored PackBits image decoder chunk. |
| [app/shiny_app/www/avivator/assets/pako.esm-Cram60i4.js](../app/shiny_app/www/avivator/assets/pako.esm-Cram60i4.js) | Vendored compression/decompression support chunk. |
| [app/shiny_app/www/avivator/assets/raw-BaiDF_Ic.js](../app/shiny_app/www/avivator/assets/raw-BaiDF_Ic.js) | Vendored uncompressed image decoder chunk. |
| [app/shiny_app/www/avivator/assets/webimage-wutFjnm1.js](../app/shiny_app/www/avivator/assets/webimage-wutFjnm1.js) | Vendored browser image-decoding support chunk. |
| [app/shiny_app/www/avivator/bootstrap.mjs](../app/shiny_app/www/avivator/bootstrap.mjs) | Authored wrapper: validates image URL before vendor import, constrains URL form behavior, and adds accessibility/XML compatibility handling. |
| [app/shiny_app/www/avivator/index.html](../app/shiny_app/www/avivator/index.html) | Viewer entry document and Content Security Policy; loads the authored bootstrap. |
| [app/shiny_app/www/avivator/manifest.json](../app/shiny_app/www/avivator/manifest.json) | Browser web-app metadata; unrelated to the scientific release manifest.json. |
| [app/shiny_app/www/avivator/robots.txt](../app/shiny_app/www/avivator/robots.txt) | Crawler instructions; not access control. |
| [app/shiny_app/www/avivator/url-policy.mjs](../app/shiny_app/www/avivator/url-policy.mjs) | Same-origin local TIFF URL validation and XML element-prefix compatibility helper. |
| [app/shiny_app/www/favicon-16.png](../app/shiny_app/www/favicon-16.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon-192.png](../app/shiny_app/www/favicon-192.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon-32.png](../app/shiny_app/www/favicon-32.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon-48.png](../app/shiny_app/www/favicon-48.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon-512.png](../app/shiny_app/www/favicon-512.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon-64.png](../app/shiny_app/www/favicon-64.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon.ico](../app/shiny_app/www/favicon.ico) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/favicon.png](../app/shiny_app/www/favicon.png) | Browser/bookmark/mobile icon asset; filename variants provide different sizes or formats. |
| [app/shiny_app/www/logo.png](../app/shiny_app/www/logo.png) | Symlink to ../../docs/logo.png; the image build recreates it. |

**Production and acceptance tooling**

| File | Role |
|---|---|
| [deploy/BASE_HARDENING.md](../deploy/BASE_HARDENING.md) | Why and how the vulnerable base server runtime was rebuilt, with validation/evidence limits. |
| [deploy/BROWSER_SAST.md](../deploy/BROWSER_SAST.md) | Scope and execution of authored-browser static analysis, with vendor coverage limitations. |
| [deploy/Containerfile](../deploy/Containerfile) | Multi-stage data-free OCI build: server hardening, R/Python dependencies, code, runtime identity, and startup defaults. |
| [deploy/EXTERNAL_ACCEPTANCE.md](../deploy/EXTERNAL_ACCEPTANCE.md) | Procedure for external HTTPS, whole-site Basic Auth, WebSocket, image byte-range, and sidecar checks. |
| [deploy/README.md](../deploy/README.md) | Authoritative production runbook: immutable artifacts, transfer, activation, proxy acceptance, recovery, and evidence. |
| [deploy/Renviron.pubapps.example](../deploy/Renviron.pubapps.example) | Nonsecret production environment reference; the Quadlet normally supplies these settings directly. |
| [deploy/activate-release.sh](../deploy/activate-release.sh) | Full preflight, UID/mount/schema/image checks, Quadlet rendering, service replacement, health verification, and automatic recovery. |
| [deploy/base-image.lock](../deploy/base-image.lock) | Digest-pinned Rocker/R base-image input; not the final app image identity. |
| [deploy/browser/README.md](../deploy/browser/README.md) | Browser acceptance setup, tested behavior, sanitized outputs, and remaining manual review. |
| [deploy/browser/check-app.mjs](../deploy/browser/check-app.mjs) | Puppeteer/axe harness for four tabs, downloads, real viewer canvas, range delivery, URL containment, and automated accessibility. |
| [deploy/browser/package-lock.json](../deploy/browser/package-lock.json) | Resolved browser-test dependency graph; belongs to test tooling, not the vendored viewer source tree. |
| [deploy/browser/package.json](../deploy/browser/package.json) | Direct Node package requirements and tool configuration for browser acceptance. |
| [deploy/build-image.sh](../deploy/build-image.sh) | Requires clean source and pinned base; builds/probes/scans an image and records immutable identity; optional registry push is a separate supported branch. |
| [deploy/container-entrypoint.sh](../deploy/container-entrypoint.sh) | Validates production feature/credential/UID/path contracts and writes approved worker environment profile before starting the server. |
| [deploy/containerignore](../deploy/containerignore) | Build-context whitelist and exclusions for scientific assets, credentials, local history, and caches. |
| [deploy/correct-syft-vim-cpe.py](../deploy/correct-syft-vim-cpe.py) | Narrow, recorded correction for the CRAN VIM/editor CPE collision; preserves unrelated components/findings. |
| [deploy/export-image.sh](../deploy/export-image.sh) | Exports the resolved immutable image as OCI with metadata and transfer checksums. |
| [deploy/generate-release-manifest.sh](../deploy/generate-release-manifest.sh) | Assembles release.env and outer checksum/evidence metadata; rejects unsafe or inconsistent release content. |
| [deploy/harden-base.sh](../deploy/harden-base.sh) | Builds the checksum-pinned server/runtime and executes upstream/native checks used by the Containerfile build stage. |
| [deploy/install-ci-duckdb.R](../deploy/install-ci-duckdb.R) | Installs the verified DuckDB version in CI's R test environment. |
| [deploy/install-production-r-deps.R](../deploy/install-production-r-deps.R) | Installs app R packages and enforces the pinned destiny/DuckDB production requirements with source integrity checks. |
| [deploy/islet-explorer.container.in](../deploy/islet-explorer.container.in) | Template for rootless user Quadlet: image/mounts/UID mapping, port, features, limits, healthcheck, and restart policy. |
| [deploy/islet-explorer.def](../deploy/islet-explorer.def) | Historical Apptainer/Singularity recipe; explicitly unsupported for reviewer production. |
| [deploy/lint-runtime-r.R](../deploy/lint-runtime-r.R) | Runtime R lint/correctness/security-pattern gate. |
| [deploy/load-image.sh](../deploy/load-image.sh) | Verifies OCI archive/metadata checksums and imports/verifies local image identity. |
| [deploy/r-worker.py](../deploy/r-worker.py) | Restores allowlisted configuration after Shiny Server resets the R environment; rejects invalid profile ownership/mode and disables user R startup files. |
| [deploy/requirements-python.lock](../deploy/requirements-python.lock) | Exact Python package version set for the runtime trajectory environment. |
| [deploy/rollback-release.sh](../deploy/rollback-release.sh) | Selects the recorded previous release and invokes the normal verified activation path. |
| [deploy/run.sh](../deploy/run.sh) | Manual Podman smoke launcher, not persistent production supervision; a pre-existing working change adds runtime health configuration. |
| [deploy/scan-browser-code.sh](../deploy/scan-browser-code.sh) | Extracts/checks authored and R-embedded browser code with the configured static-analysis toolchain. |
| [deploy/scan-image.sh](../deploy/scan-image.sh) | Generates SBOM/dependency/security evidence for an exact image and enforces the project High/Critical gate. |
| [deploy/shiny-server.conf](../deploy/shiny-server.conf) | Shiny prefix and explicit static TIFF routes, error sanitization, worker logs, timeouts, and disabled directory listings. |
| [deploy/tests/smoke-shiny-server.py](../deploy/tests/smoke-shiny-server.py) | Real Shiny Server/R/static-route/environment smoke used during hardened server construction. |
| [deploy/tests/test-deployment-recovery.sh](../deploy/tests/test-deployment-recovery.sh) | Exercises activation failure and previous-service recovery contracts with controlled fixtures. |
| [deploy/tests/test-deployment-verifier.sh](../deploy/tests/test-deployment-verifier.sh) | Exercises local deployment verifier failure/success conditions. |
| [deploy/tests/test-image-evidence.sh](../deploy/tests/test-image-evidence.sh) | Exercises image/release evidence and identity contracts. |
| [deploy/tests/test-viewer-boundary.mjs](../deploy/tests/test-viewer-boundary.mjs) | Checks URL containment, missing/external/traversal inputs, XML compatibility, CSP order, and vendor hashes. |
| [deploy/tests/test_reviewer_access.py](../deploy/tests/test_reviewer_access.py) | Synthetic transport tests for external acceptance logic; not proof of the real UF proxy. |
| [deploy/tests/test_sbom_correction.py](../deploy/tests/test_sbom_correction.py) | Regression tests bounding the single CRAN VIM/editor CPE correction. |
| [deploy/verify-deployment.sh](../deploy/verify-deployment.sh) | Verifies local health, published-port boundary, UID/mounts, feature flags, dependencies, immutable image, and optional Range behavior. |
| [deploy/verify-release.sh](../deploy/verify-release.sh) | Verifies assembled release hashes, image count, metadata, and verification stamp; full and quick modes differ. |
| [deploy/verify-reviewer-access.py](../deploy/verify-reviewer-access.py) | Checks actual external HTTPS/auth/WebSocket/sidecar/Range boundary with credential-file handling and sanitized evidence. |

**Operating and research documentation**

| File | Role |
|---|---|
| [docs/HANDOFF.md](../docs/HANDOFF.md) | Historical development and statistical-correctness handoff; old deployment instructions are superseded. |
| [docs/admin_guide.md](../docs/admin_guide.md) | Data rebuild and administration reference; defers production authority to deploy/README and the UFRC handoff. |
| [docs/competitive_landscape.md](../docs/competitive_landscape.md) | Research/product comparison background; not deployment configuration. |
| [docs/pubapps_architecture.png](../docs/pubapps_architecture.png) | Rendered request architecture diagram. |
| [docs/pubapps_architecture.svg](../docs/pubapps_architecture.svg) | Editable source for the request architecture diagram. |
| [docs/pubapps_log_inventory.md](../docs/pubapps_log_inventory.md) | Log streams, minimization rules, collection/retention questions, and responsibilities. |
| [docs/reviewer_guide.md](../docs/reviewer_guide.md) | Instructions for the four-tab reviewer profile, access, exports, and limitations. |
| [docs/ufrc_readiness.md](../docs/ufrc_readiness.md) | Requirement-to-evidence index, recorded release results, unresolved limitations, and remaining external acceptance. |
| [docs/user_guide.md](../docs/user_guide.md) | Detailed application feature reference, including retained optional/development features. |

**Scientific preparation scripts**

| File | Role |
|---|---|
| [scripts/CLAUDE.md](../scripts/CLAUDE.md) | Scientific preparation scripts, execution order, known invariants, and development history. |
| [scripts/Comb_Xen_Phen.ipynb](../scripts/Comb_Xen_Phen.ipynb) | Exploratory Xenium/PhenoCycler combination notebook; upstream research, not reviewer runtime. |
| [scripts/annotate_ome_tiff_channels.py](../scripts/annotate_ome_tiff_channels.py) | Preparation tool that changes OME-XML channel names/colors in-place; use on source/staging, not a frozen release. |
| [scripts/auto_color_ome_tiff_channels.py](../scripts/auto_color_ome_tiff_channels.py) | Preparation tool that assigns channel colors in TIFF metadata; can change files in-place. |
| [scripts/build_barlow_features.py](../scripts/build_barlow_features.py) | Builds core/peri-islet composition features consumed by runtime diffusion/PAGA analysis. |
| [scripts/build_donor_metadata.py](../scripts/build_donor_metadata.py) | Builds curated donor metadata for downstream application/data preparation. |
| [scripts/build_h5ad_for_app.py](../scripts/build_h5ad_for_app.py) | Assembles the development app H5AD from scientific, annotation, donor, neighborhood, and clustering inputs. |
| [scripts/build_master_excel.py](../scripts/build_master_excel.py) | Builds the legacy/development summary workbook from upstream exported tables. |
| [scripts/build_release_data.R](../scripts/build_release_data.R) | Creates the immutable runtime bundle: prepared Parquet summaries, partitioned tissue, indexed geometry, selected copied files, images, and manifest. |
| [scripts/clear_shiny_cache.sh](../scripts/clear_shiny_cache.sh) | Legacy local cache/restart helper; not the current release activation procedure. |
| [scripts/compile_donors.py](../scripts/compile_donors.py) | Compiles donor-level source metadata/tables for upstream preparation. |
| [scripts/compute_neighborhood_metrics.py](../scripts/compute_neighborhood_metrics.py) | Computes peri-islet composition, enrichment, immune/spatial metrics, and distance tables from single-cell inputs. |
| [scripts/compute_pseudotime.R](../scripts/compute_pseudotime.R) | Older workbook-based pseudotime analysis script; not the active runtime PAGA entrypoint. |
| [scripts/compute_split_pseudotime.py](../scripts/compute_split_pseudotime.py) | Historical count-stratified pseudotime preparation that writes H5AD data; active UI now uses the documented composition/count-range model. |
| [scripts/convert_to_ometiff.sh](../scripts/convert_to_ometiff.sh) | Upstream TIFF-to-pyramidal-OME-TIFF conversion wrapper; optional Docker/JAR paths are preparation alternatives, not PubApps runtime choices. |
| [scripts/convert_to_omezarr.sh](../scripts/convert_to_omezarr.sh) | Upstream OME-Zarr conversion helper; selected reviewer image format is OME-TIFF. |
| [scripts/derive_gates.py](../scripts/derive_gates.py) | Derives/calibrates per-image marker gates and records scientific provenance. |
| [scripts/diagnose_islets.R](../scripts/diagnose_islets.R) | Inspects workbook/islet preparation and count/key behavior. |
| [scripts/diagnose_marker_regions.R](../scripts/diagnose_marker_regions.R) | Inspects marker-region tables and data-preparation assumptions. |
| [scripts/extract_browser_scripts.R](../scripts/extract_browser_scripts.R) | Parses literal browser programs embedded in R without evaluating the R source; supplies SAST inputs. |
| [scripts/extract_per_donor_tissue.py](../scripts/extract_per_donor_tissue.py) | Builds per-donor all-cell tissue CSVs for later release conversion/spatial lookup. |
| [scripts/extract_per_islet_cells.py](../scripts/extract_per_islet_cells.py) | Builds per-islet cell CSVs for cell-level drilldown; old header count/filter examples are historical. |
| [scripts/fix_tcell_subtypes.py](../scripts/fix_tcell_subtypes.py) | Scientific correction of CD4/CD8 subtype assignments in canonical single-cell data; changes scientific inputs. |
| [scripts/get_bioformats2raw.sh](../scripts/get_bioformats2raw.sh) | Downloads upstream conversion JARs for preparation workflows. |
| [scripts/install_avivator.sh](../scripts/install_avivator.sh) | Legacy viewer installation/mirroring helper with an old shiny/www path; not the current reviewed/pinned vendor update procedure. |
| [scripts/install_shiny_deps.R](../scripts/install_shiny_deps.R) | Main R dependency installer/manifest used by development and wrapped by the production installer. |
| [scripts/inventory_release_data.py](../scripts/inventory_release_data.py) | Produces sanitized schema/count/asset inventory without serializing donor rows or embedded text; supports disclosure review but cannot establish classification alone. |
| [scripts/make_offsets_json.py](../scripts/make_offsets_json.py) | Creates TIFF offset sidecars for efficient browser access. |
| [scripts/na_audit.R](../scripts/na_audit.R) | Legacy/development missing-value audit for workbook tables. |
| [scripts/plot_scvi_comparison.py](../scripts/plot_scvi_comparison.py) | Generates supplemental comparison figures for original and retrained single-cell models. |
| [scripts/reaggregate_islets.py](../scripts/reaggregate_islets.py) | Reaggregates single cells to islets and computes scientific trajectory/clustering outputs upstream of the app release. |
| [scripts/recompute_viz_umap.py](../scripts/recompute_viz_umap.py) | Recomputes a visualization embedding and rebuilds affected app scientific inputs. |
| [scripts/retrain_scvi_no_covariates.py](../scripts/retrain_scvi_no_covariates.py) | Research comparison training without age/gender covariates; not a reviewer runtime requirement. |
| [scripts/run_ci_r_tests.R](../scripts/run_ci_r_tests.R) | Runs selected R suites, retains machine-readable outcomes, and rejects skipped/failing CI checks. |
| [scripts/scvi_comparison_results/FigureGuide.txt](../scripts/scvi_comparison_results/FigureGuide.txt) | Interpretive guide to the supplemental model-comparison figure. |
| [scripts/scvi_comparison_results/comparison_results.json](../scripts/scvi_comparison_results/comparison_results.json) | Recorded model-comparison metrics. |
| [scripts/scvi_comparison_results/diagnostics_log.txt](../scripts/scvi_comparison_results/diagnostics_log.txt) | Recorded scientific model-comparison diagnostic output. |
| [scripts/scvi_comparison_results/islet_metadata.csv](../scripts/scvi_comparison_results/islet_metadata.csv) | Islet metadata aligned with the comparison embeddings. |
| [scripts/scvi_comparison_results/new_islet_latent.npy](../scripts/scvi_comparison_results/new_islet_latent.npy) | Numeric array of comparison-model islet latent values. |
| [scripts/scvi_comparison_results/orig_islet_latent.npy](../scripts/scvi_comparison_results/orig_islet_latent.npy) | Numeric array of original-model islet latent values. |
| [scripts/scvi_comparison_results/supplemental_scvi_comparison.png](../scripts/scvi_comparison_results/supplemental_scvi_comparison.png) | Rendered supplemental model-comparison figure. |
| [scripts/test_shiny_prep.R](../scripts/test_shiny_prep.R) | Legacy/development scientific preparation smoke/audit script. |
| [scripts/tests/test_browser_extraction.R](../scripts/tests/test_browser_extraction.R) | Tests safe extraction of browser literals without executing arbitrary R. |
| [scripts/tests/test_inventory_release_data.py](../scripts/tests/test_inventory_release_data.py) | Tests inventory behavior/minimization with controlled fixtures. |
| [scripts/tests/test_release_parity.R](../scripts/tests/test_release_parity.R) | Tests data-conversion/parity edge cases, including blank strings and numeric missingness. |
| [scripts/verify_release_data.R](../scripts/verify_release_data.R) | Checks data manifest/schema/artifacts/hashes/row counts before release promotion; core/full modes serve different stages. |
| [scripts/verify_release_parity.R](../scripts/verify_release_parity.R) | Compares real canonical source preparation and converted release for scientific parity with sanitized evidence. |
| [scripts/verify_static_images.sh](../scripts/verify_static_images.sh) | Legacy local Shiny/static-route diagnostic; use deploy/verify-deployment and external acceptance for the current production contract. |

**CODEX research pipeline adapters**

| File | Role |
|---|---|
| [scripts/codex_phenocycler/build_cells_parquet.py](../scripts/codex_phenocycler/build_cells_parquet.py) | CODEX-specific QuPath TSV ingestion/partitioned Parquet preparation; adapts upstream cohort assumptions. |
| [scripts/codex_phenocycler/codex_floors.py](../scripts/codex_phenocycler/codex_floors.py) | CODEX-specific marker floor estimation used in the research port. |
| [scripts/codex_phenocycler/codex_lineage.py](../scripts/codex_phenocycler/codex_lineage.py) | CODEX-panel broad lineage configuration/logic. |
| [scripts/codex_phenocycler/codex_panel.py](../scripts/codex_phenocycler/codex_panel.py) | CODEX marker/channel/cohort compatibility definitions and checks. |
| [scripts/codex_phenocycler/codex_positivity.py](../scripts/codex_phenocycler/codex_positivity.py) | Converts reviewed normalization references into per-cell positivity and compartment comparisons. |
| [scripts/codex_phenocycler/codex_redsea.py](../scripts/codex_phenocycler/codex_redsea.py) | CODEX adapter for upstream REDSEA lateral spillover correction. |
| [scripts/codex_phenocycler/codex_reference_selection.py](../scripts/codex_phenocycler/codex_reference_selection.py) | Builds target-reference review dossiers; not a production phenotype writer. |
| [scripts/codex_phenocycler/codex_restore_pair_review.py](../scripts/codex_phenocycler/codex_restore_pair_review.py) | Builds CODEX-specific RESTORE pair image-review bundles through the upstream implementation. |
| [scripts/codex_phenocycler/codex_restore_validation.py](../scripts/codex_phenocycler/codex_restore_validation.py) | CODEX adapter for upstream RESTORE validation/scoring. |
| [scripts/codex_phenocycler/docs/HANDOFF.md](../scripts/codex_phenocycler/docs/HANDOFF.md) | Research-port state, decisions, and unresolved scientific review work; separate from PubApps deployment. |
| [scripts/codex_phenocycler/docs/codex_restore_pairs.md](../scripts/codex_phenocycler/docs/codex_restore_pairs.md) | Scientific pair/reference definitions and their rationale. |
| [scripts/codex_phenocycler/docs/rejected/adaptive_restore_divisor.patch](../scripts/codex_phenocycler/docs/rejected/adaptive_restore_divisor.patch) | Rejected research implementation retained for provenance; not applied production code. |
| [scripts/codex_phenocycler/export_cells_geojson.groovy](../scripts/codex_phenocycler/export_cells_geojson.groovy) | QuPath cell polygon export for upstream scientific preprocessing. |
| [scripts/codex_phenocycler/export_cells_geojson.sh](../scripts/codex_phenocycler/export_cells_geojson.sh) | Wrapper for the QuPath cell-geometry export workflow. |
| [scripts/codex_phenocycler/test_codex_panel.py](../scripts/codex_phenocycler/test_codex_panel.py) | Checks CODEX marker/channel/donor mapping against relevant scientific inputs. |
| [scripts/codex_phenocycler/test_codex_restore_validation.py](../scripts/codex_phenocycler/test_codex_restore_validation.py) | Checks the CODEX adapter's configuration and patching contracts. |
| [scripts/codex_phenocycler/verify_compartment_means.py](../scripts/codex_phenocycler/verify_compartment_means.py) | Checks QuPath compartment exports are complete and aligned before downstream normalization/review. |
| [scripts/codex_phenocycler/verify_geojson.py](../scripts/codex_phenocycler/verify_geojson.py) | Checks polygon IDs, coordinates, cohort, and cell/image correspondence before spillover processing. |

**Data provenance, research, and external pipeline**

| File | Role |
|---|---|
| [data/DATA_PROVENANCE.md](../data/DATA_PROVENANCE.md) | Canonical scientific lineage and historical validation, with a current prepared-table count correction. |
| [data/donor_metadata.csv](../data/donor_metadata.csv) | Curated donor clinical/demographic metadata consumed by the application and release builder; Git-tracked exception to bulk-data ignores. |
| [data/phenotype_gates_per_image.csv](../data/phenotype_gates_per_image.csv) | Per-image marker gate values used in scientific preparation/provenance. |
| [data/phenotype_gates_provenance.csv](../data/phenotype_gates_provenance.csv) | Provenance for the per-image marker gates. |
| [external/Phenocycler_Analysis](../external/Phenocycler_Analysis) | Git submodule reference to an upstream research pipeline; source tree is separately versioned and not part of the Shiny image. |
| [islet_analysis/CLAUDE.md](../islet_analysis/CLAUDE.md) | Islet aggregation/research pipeline notes; the only top-level Git-tracked entry within that otherwise ignored analysis directory. |
| [notebooks/rebuild_trajectory.ipynb](../notebooks/rebuild_trajectory.ipynb) | Scientific trajectory rebuild/validation notebook; historical counts must be interpreted in context. |
| [notebooks/scvi_qc_validation.ipynb](../notebooks/scvi_qc_validation.ipynb) | Single-cell latent-space and donor/cell-type QC notebook. |
| [single_cell_analysis/CODEX_Panc_Analysis.ipynb](../single_cell_analysis/CODEX_Panc_Analysis.ipynb) | Upstream pancreatic CODEX single-cell analysis notebook. |
| [single_cell_analysis/CODEX_Panc_scimap_Analysis.ipynb](../single_cell_analysis/CODEX_Panc_scimap_Analysis.ipynb) | Upstream single-cell phenotyping/spatial analysis notebook using scimap. |
| [single_cell_analysis/SCVI_CODEX_copy.ipynb](../single_cell_analysis/SCVI_CODEX_copy.ipynb) | Alternate/historical scVI analysis notebook copy; not runtime. |
| [single_cell_analysis/SCVI_CODEX_v2.ipynb](../single_cell_analysis/SCVI_CODEX_v2.ipynb) | Upstream scVI single-cell analysis/training notebook; not runtime. |
| [single_cell_analysis/phenotype_rules1.csv](../single_cell_analysis/phenotype_rules1.csv) | Earlier single-cell phenotyping rule reference. |

**Historical files**

| File | Role |
|---|---|
| [archive/README.md](../archive/README.md) | Explains why archived code/notes are retained and which files remain active elsewhere. |
| [archive/docs/HANDOFF_2026-07-25_shiny-correctness-pass.md](../archive/docs/HANDOFF_2026-07-25_shiny-correctness-pass.md) | Historical correctness-pass record, superseded by later handoffs. |
| [archive/misc/FileReview.txt](../archive/misc/FileReview.txt) | Historical implementation summary. |
| [archive/misc/shiny-Ex.R](../archive/misc/shiny-Ex.R) | Old generated Shiny package example, not application runtime code. |
| [archive/phenotype_rules2.csv](../archive/phenotype_rules2.csv) | Superseded phenotyping rule set. |
| [archive/shiny_app/app_modular.R](../archive/shiny_app/app_modular.R) | Historical snapshot of the modular entrypoint; not loaded. |
| [archive/shiny_app/app_original.R](../archive/shiny_app/app_original.R) | Historical pre-modular monolithic app; not loaded. |

**Ignored/local families and what their individual files mean**

| Family | Role and relationship to production |
|---|---|
| `data/cells/<case>_Islet_<number>.csv` | One small cell-level drilldown table per islet; selected release copies these. |
| `data/donors/<case>.csv` | Donor tissue source table; release builder converts it to partitioned Parquet. |
| `data/json/*.geojson`, `data/gson/*.geojson.gz` | Segmentation sources in plain/compressed representations; builder creates indexed GeoPackages. |
| `data/gson/gson_cache/` | Local cached geometry/preprocessing; not a substitute for generated release geometry. |
| `data/parquet/` | Local prepared/experimental Parquet layout; presence alone does not make a manifest-valid production release. |
| `data/islet_explorer.h5ad`, `master_results.xlsx` | Development/canonical preparation inputs; the production summary loader uses the generated release. |
| `data/adata_ins_root.h5ad` and selected lookup/feature/rule CSVs | Scientific runtime inputs copied by the release builder. |
| `data/annotations.tsv`, neighborhood/distance CSVs | Source/intermediate scientific tables; builder selects/transforms the required subset. |
| `*.ome.tiff` in `www/local_images/` | Full-resolution/pyramidal tissue image, one per donor in the selected cohort. |
| `*.ome.xml` | Image/channel metadata sidecar; included in the authorized image disclosure and authenticated routes. |
| `*.ome.tiff.offsets.json`, `*.offsets.json` | TIFF offset indexes and compatibility symlinks for browser random access. |
| `.Renviron`, `.Rprofile`, `.Rhistory`, `.RData` | Local environment/startup/history/workspace state; not the reviewer runtime configuration authority. The production worker disables user startup files. |
| `APP_LOADED_FLAG` files | Local historical runtime markers; excluded from production build context. |
| Root donor workbook and `CellCompartments.tsv` | Source clinical metadata and large QuPath export; not a reason to transfer the entire repository to PubApps. |
| `islet_analysis/` H5AD/TSV/CSV/PNG/notebooks | Research aggregations, validation, figures, and intermediate artifacts. |
| `single_cell_analysis/` H5AD/TSV | Upstream single-cell data/model inputs and products; not runtime training on PubApps. |
| Executed notebooks and generated figures | Research execution records/outputs, distinct from versioned notebook source. |
| `scripts/Python/`, `scripts/R_files/`, `scripts/Shell_scripts/` | Ignored organizational copies and legacy experiments. Several have stale paths or mutable deployment behavior; individually review before reuse. |
| `external/Phenocycler_Analysis/phenocycler/` | Upstream reusable research package: ingestion, correction, lineage, validation, figures, parallelism, and configuration. |
| External `notebooks/`, `scripts/groovy/`, `scripts/slurm/`, `scripts/R/`, `tests/` | Upstream research entrypoints, QuPath jobs, scheduled processing, diagnostics, and tests. Not Shiny production services. |
| Nested `external/RESTORE/` | Removed 2026-10-15: deinitialized and excluded from recursive updates (`submodule.external/RESTORE.update none`). Nothing in this repository used it; the empty gitlink directory remains only because the pinned `8db37af` upstream still declares it. See [archive/README.md](../archive/README.md). |
| `tmp/ufrc-release-20260908/artifacts/` | Multiple retained candidates, including rejected/historical data and the assembled selected release; identify by manifest, not directory size. |
| `tmp/ufrc-release-20260908/evidence/` | Local restricted build/scan/parity/browser/request evidence. Durable retention must be arranged outside disposable local staging. |
| `node_modules/`, `.venv/`, `__pycache__/` | Installed third-party dependencies and generated interpreter caches. Inventory by environment/manifest and path; these are not individually authored app modules. |
| `literature/`, `*.pdf` | Local scientific references and outputs. |
| `.vscode/`, local assistant settings | Workstation preferences/tooling; not reviewer configuration. |
| `.git/`, submodule `.git` reference files | Version-control object/index/reference metadata. Source-file review does not imply rewriting Git history. |
| `.agents/`, `.codex/` | Managed session/tooling directories; empty at the workspace root during inspection. |

Binary scientific assets and generated release evidence cannot usefully undergo a prose cleanup. Review them for their intended role and metadata; preserve approved frozen artifacts. Any later changes to runtime files or hash-bound viewer assets must become a new verified candidate.

**Suggested review batches**

Start with the runtime entrypoint and shared helpers, then the four visible module pairs, browser assets, data/release builders, container/startup/service files, security/acceptance tooling, documentation, and finally research/archive/tooling material. The review register starts every entry as `pending_user_review`; no AI-related deletion has been applied.
