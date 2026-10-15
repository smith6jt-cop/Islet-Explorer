# Commit preparation — 2026-09-17

Status: commit plan prepared; no commit or push performed. The real Git index
remains unchanged. A temporary index is used to inspect the proposed source-only
patch and preserve any subsequent user staging choices.

## Proposed commit

`Add 21-donor candidate workflow and trajectory review safeguards`

Include the pending implementation as one coherent development checkpoint:

- Explicit cohort/source selection, corrected donor preparation, immutable
  candidate builders, native QuPath review export and audit scripts.
- Fixed offline trajectory, complete 31-marker heatmap, donor-filtered outputs,
  validated ratio/count contracts, islet-identity review requirements and tests.
- Offline hypothesis runners/report packaging, removal of research from the
  current app, schema-3 manifests and schema-1/2 compatibility.
- Preview/deployment/CI changes, current provenance and readiness records,
  dated implementation notes, and archived source for the former trajectory.
- Updated `.gitignore` and the required static `data/phenotype_rules.csv` rule
  definition. The existing tracked donor/gate provenance tables are unchanged.

This records implemented code and reproducible diagnostics. It does **not**
claim accepted islet segmentation, validated phenotype gates, available real
trajectory ratios, accepted research findings or UFRC deployment acceptance.
The owner's HPAP-078 finding remains an explicit scientific defect to resolve.

## Ignore audit

The existing `data/*` rule already excludes all candidate datasets, native
review projects, generated charts/reports and evidence. Added patterns also
exclude scientific matrix/image/geometry files outside that tree, native QuPath
projects/object files and their backups, root-level generated reports, runtime
logs/PIDs, Python/R caches and local credential files. `.env.example` and
`.env.template` remain eligible as templates.

`data/phenotype_rules.csv` was previously ignored despite being a required
builder/test input. It is a 1,322-byte static rule table and is now explicitly
allowed. New cohort CSVs and the hypothesis registry remain eligible source.
The new source cohort CSV's line endings were normalized to LF for Git's
whitespace check. All 21 parsed donor rows and the candidate's copied cohort
bytes remain unchanged; the first preflight diagnostic is preserved separately.
No images, cell tables, object hierarchies or numerical result artifacts should
appear in the proposed patch. Ignore rules do not untrack files already in Git,
so the preflight also checks tracked changes and the candidate index.

## Review and commit sequence

1. Inspect the exact proposed path list, patch statistics and patch under the
   ignored `data/development/all-donor-20260917/evidence/commit-preflight/` folder.
   The preflight records source hashes, checks ignored/allowed examples, scans
   pending text for common credential signatures without printing values, and
   reports file sizes. This is a focused scan, not proof that every possible
   secret format has been detected.
2. Confirm all changed files belong to this checkpoint, including earlier
   pending cohort/deployment documentation. Do not stage unrelated later edits.
   Use the reviewed explicit path list rather than a blanket `git add .`.
3. Verify the recorded R tests/lint, Python pipeline tests, synthetic deployment
   checks and current browser evidence. Scientific input review remains pending;
   test success must not change that status.
4. Stage only the reviewed paths, inspect `git diff --cached --check` and
   `git diff --cached --stat`, then create the proposed commit. A source-hash
   comparison must be repeated if files changed after the preflight snapshot.
5. Record the resulting commit ID in a new evidence record, keeping the existing
   source snapshot and frozen release evidence intact. Pushing and deployment
   are separate actions.

The current app and scientific limitations are documented in the
[revision record](trajectory_revision_20260917.md). The restricted preflight
files and patch are local review material and are intentionally excluded from
the proposed commit.
