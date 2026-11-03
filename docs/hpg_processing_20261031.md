# Typing processing on HiPerGator (owner decisions of 2026-10-31)

**Subsequent selected-component build, November 2:** the completed v2 app and
QuPath builds use scientific source `67ee6c9`; the matched image/release builder
uses `e5689f6`. Production `pubapps-2026-11-components-v1` is now live and
verified (2026-11-03 02:42 UTC); dev remains stopped and masked. Read [the v2 build record](component_build_v2_20261102.md),
[resolution contract](component_resolution_20261102.md), and
[actual deployment status](container_security_20261102.md) for current work.
The version-specific source paths, failed diagnostics and prior service state
below remain historical evidence; they do not reopen the delegated selection
or production decision.

**2026-11-02 correction:** the old `run_v20.py` full-calibration phase is historical
and must not be used for a Lymphatic-only revision. Use the explicit calibration/context
baseline and automatic invariant audit in [the workflow guide](workflow_user_guide.md).
The [component comparison](v19_v20_component_comparison_20261102.md) records the real-data
regression and corrected diagnostic replay. Local replay evidence does not replace a fresh
authoritative HPG candidate. Source synchronization uses the explicit changed-file rule below;
sealed scientific inputs and the four protected modules are never overwritten.

The remaining steps of [the 2026-10-31 typing handoff](typing_handoff_20261031.md) run on HiPerGator (HPG) from a
copy of the workstation's data. Owner decisions, 2026-10-31:

- **Scope.** The whole pipeline moves to HPG, not only detection.
- **Guards.** On HPG a guard accepts a copied file while the file still matches a seal. The seal is recorded after
  the copy, from the file's sha256, size and modification time.
- **QOS.** The burst QOS `maigan-b` may be used for this work, an exception to the 2026-10-01 rule (investment QOS
  only).

Approved plan: `/home/smith6jt/.claude/plans/review-the-handoff-and-snazzy-dijkstra.md`.

## Why a seal

The typing guards record each protected file's `{path, size, mtime_ns, ctime_ns, inode}` and compare it exactly:
- `structure_generation._stat`;
- `new_detection_filter.file_stat`;
- `refresh_ready_cells.file_stat`;
- `run_remaining_qc_pipeline.stat_identity`;
- `final_reconciliation`;
- `collagen_structures`' raw-image size and mtime.

No copy can carry those values over:
- Globus and rsync keep only the modification time.
- Setting a timestamp stamps a new ctime (`man 2 utimensat`).
- Each filesystem assigns its own inodes.
- HPG's Lustre `/blue` keeps mtime to the whole second.

The guard functions themselves can't be edited either. Four of them sit in files whose own sha256 is pinned in
records:
- `new_detection_filter.py`, in 20 filter manifests;
- `refresh_ready_cells.py`, in 80 REDSEA cache records;
- `collagen_structures.py`;
- `verify_display_source_pixels.py`.

Instead, `scripts/hpg/sitecustomize.py` (HPG only, on `PYTHONPATH` inside the site container) wraps `os.stat`
(`scripts/hpg/site_identity.py`):
- **Translated.** A regular file whose current HPG (inode, size, mtime_ns, ctime_ns) equals a sealed tuple is
  reported with its workstation inode, device and times.
- **Everything else is reported as is.** That covers a changed, replaced or new file, a directory, `lstat`,
  `fstat`, `DirEntry.stat` and `dir_fd` calls. A change after sealing therefore fails the guards exactly as on the
  workstation.
- **A process without the hook fails the guards** (`python -S`, `-I`, `-E`, or an unset `CODEX_SITE_SEAL`).
- **The hook is pinned.** It refuses a seal whose sha256 differs from `CODEX_SITE_SEAL_SHA256` and ends the
  process with exit code 78.
- **Journals record it.** `run_v20.py` and `activate_pdpn_gate_selection.py` record it in every journal
  (`site_identity`), and refuse to run under SLURM without it.
- **Pre-edit copies are kept**: `*.before_site_identity.py`, which match the sha256 the 10-30 journals and
  `activation.json` record.

## Layout

| HPG | Seen inside the container as |
|---|---|
| `/blue/maigan/smith6jt/islet-explorer-hpg/root` | `/data/panc_CODEX/Islet-Explorer` (one bind mount: activation renames across `$RUN` and `output/`) |
| `/blue/maigan/smith6jt/islet-explorer-hpg/qpi` | `/data/QuPath_Projects/QuPath_project_IsletExplorer` (read-only) |
| `/home/smith6jt/QuPath` | the same (QuPath 0.7.0, core jars sha256-identical to the workstation's) |
| `/home/smith6jt/miniconda3/envs/scvi-env` | a symlink to `/blue/maigan/smith6jt/miniforge3/envs/scvi-env` (same package versions) |

- `scripts/hpg/run_in_site.sh` runs a command in the container: `apptainer exec --cleanenv`, with the host root as
  the container and the binds above.
- `scripts/hpg/jobs/site_job.sbatch` wraps it for SLURM on `hpg-turin`. Those nodes are Zen 5 with AVX-512, the
  same instruction-set family as the workstation's Zen 4.
- QuPath's OpenCV natives need the `lib_stdcxx` symlink (`hpg_env.sh`).

## Rules

1. **Capture once.** `site_seal.py capture` runs on the workstation with QuPath closed and
   `run_remaining_qc_pipeline.py --phase ready` passing. After it, the data in the workstation's `$RUN`,
   `output/typing_20261030` and `output/typing_resume_20261029` is frozen, and only HPG writes data.
2. **Code is edited on the workstation only.** Commit it and deliver it to HPG as an explicit changed-file list plus
   `.git`.
   - No whole-tree rsync into the sealed tree.
   - No `git checkout`, `stash` or `reset` there.
   - `redsea_cells.py`, `codex_typing.py`, `refresh_ready_cells.py` and `refresh_filtered_cells.py` have pinned
     inodes and must never be rewritten on either side.
3. **Results come back once.** They return to the workstation's identical absolute paths, which must not exist
   beforehand. Their `SHA256SUMS` is verified there.
4. **Reseal after admin changes.** An `lfs migrate` or an administrator's `chmod`/`chown` on `/blue` changes ctime,
   so the copy then reports its real identity and the gates fail. Run `verify-and-seal` again into a new seal and pin
   its sha256.
5. **Owner rules hold.** No interpretation. New outputs go to fresh destinations. Every gating geometry is checked
   against the image on every image before anything downstream runs.

## Commands

```bash
# workstation
python scripts/hpg/site_seal.py capture --list LIST --out output/hpg_transfer_20261031/capture --git-root .
# HPG, from a login node (ssh hpg): smoke test, seal, gates
cd /blue/maigan/smith6jt/islet-explorer-hpg/root
/opt/slurm/bin/sbatch --job-name=smoke scripts/hpg/jobs/site_job.sbatch --no-seal bash scripts/hpg/smoke.sh
/opt/slurm/bin/sbatch --job-name=seal --cpus-per-task=32 --mem=64G --time=03:00:00 scripts/hpg/jobs/site_job.sbatch \
    --no-seal python scripts/hpg/site_seal.py verify-and-seal --identities output/hpg_transfer_20261031/capture/workstation_identities.json \
    --identities-sha256 SHA --out output/hpg_transfer_20261031/hpg_site_seal.json
/opt/slurm/bin/sbatch --job-name=gates --cpus-per-task=8 --mem=64G --time=02:00:00 scripts/hpg/jobs/site_job.sbatch \
    python scripts/hpg/gate_dump.py --out output/hpg_transfer_20261031/gate_dump_hpg.json
```

The steps, their outcomes and the measured equivalence with the workstation are recorded below as they run.

## Record

### 2026-11-01

**Before the copy.**
- The despeckled Lymphatic detection D finished on the workstation on 10-31 (20 of 20).
- The evaluation and the D/DG renders (`output/typing_20261030/lymphatic_despeckle/renders_v2`, six tiles per image,
  every image viewed, notes in `review.md`) were reported. The owner chose **DG** ("Do dg").

**Workstation no-op check.** `scripts/hpg/gate_dump.py` was run before and after the hook commit `0ca1dfd`: code
inventory, REDSEA readiness (all 20 current), `structure_generation.validate` and the `run_v20` plans were
byte-identical.

**Capture and copy.**
- The capture (`output/hpg_transfer_20261031/capture`) holds 21,421 files (570.7 GB) and 57 symlinks.
  `workstation_identities.json` has sha256 434295bf46895ced7783607414b881fce8f67a3eeac15365c1f5323516a00ecd.
- The list was built by `build_transfer_list.py` from three sources: the folders the steps read, every existing file
  named in the run records, and the files the read-only gates opened (`audit_opens.py`). It excludes:
  - the nPOD spreadsheet;
  - donor 78;
  - the private blind keys;
  - `filter/.staging`;
  - the raw run's `tmp`/`work`;
  - bytecode caches;
  - other `output/` folders.
- Copied with `rsync -aH` (`transfer.sh`) at about 111 MB/s, 1 Gb/s line rate: the QuPath folder (60 files,
  163.4 GB) and the repository part (21,361 files, 407.4 GB, 57 symlinks), both exit 0.

**Container.** Apptainer refused the host root as a sandbox (smoke job 48752194). The site runs in a Rocky Linux 9
image with the host `/usr` bound read-only. Smoke job 48753532 passed on hpg-turin:
- Python 3.13.9 with the same package versions;
- git at the copied commit;
- QuPath 0.7.0 with its OpenCV natives;
- `flock` blocking a second process;
- renames across `data/` and `output/`;
- the QuPath folder read-only.

**Seal.** Job 48766607 verified all 21,421 files (sha256, size, mtime to the second) and 57 symlinks in 69 s. Seal
`output/hpg_transfer_20261031/hpg_site_seal.json`, sha256
54cad8c3c38e80eaf2405c3df34e4fac3168901b9fe892959691e56955f52ef0 (pinned in `hpg_env.sh`, commit `4afd6f9`).

**Gate checks on HPG** (on the copied code, before any code delivery):
- **With the hook** (job 48766676): REDSEA readiness for all 20, the code inventory, the `run_v20` plans, and
  `structure_generation.validate`'s binding, manifest, provenance, paths and 6,246 hash-pinned files are identical
  to the workstation dump.
- **The one difference** is the stat values of 62 directories, which the hook leaves untranslated. Their entry lists
  are identical, and directories are never compared with recorded values. The 20 source-image identities match
  through the hook.
- **Without the hook** (job 48766677): fails at `new_detection_filter.py:412`, "Aab_6450: immutable raw source
  changed".

**Reproduction on HPG** (equivalence check `scripts/hpg/equivalence.py`, job 48771468, output
`output/hpg_equivalence_20261031/equivalence.json`):
- **Despeckled PDPN:** unmodified copies of `despeckle_pdpn.py` (job 48767362) gave 20 of 20 rasters
  byte-identical to the workstation's recorded sha256.
- **Detection:** `lymphatic_despeckle_run.py` gave D for Aab_6505 and ND_6356 and original mode for Aab_6505, all
  three equal in normalized geometry. The original-mode Aab_6505 objects equal the reviewed re-run objects.
- **Phases:** `run_v20.py --phase structures` then `--phase typing` on the installed option-1 selection (job
  48767711; clive, 24 cores, 300 GB; 31 min) gave 82 of 82 structures outputs and 95 of 95 typing outputs
  byte-identical to the 10-30 workstation journals (`runs/20261030T210113.520762Z`, `runs/20261030T212756.342638Z`).
  The HPG journals are `runs/20261101T154028.006372Z` and `runs/20261101T160406.566411Z`.

**QOS.** The owner noted that the investment allocations (clive, maigan) have room, and jobs were sized to fit them
rather than using burst.

**DG install and rerun on HPG** (job 48771552; clive, 24 cores, 300 GB; 36 min):
1. `lymphatic_despeckle_gate.py` wrote `output/typing_20261030/lymphatic_despeckle_gate_v1`. All 20 gate records
   passed, and every render was byte-identical to the renders viewed before the owner's decision.
2. `activate_despeckle_gate_selection.py --activate` validated it as `despeckled_level_gate` (`active_validated`;
   6,353 selected inputs; 82 identities). Option 1's binding is kept as `selection/previous_binding.json`.
3. `run_v20.py --phase structures` (`runs/20261101T161722.560230Z`), then `--phase typing`
   (`runs/20261101T164221.481579Z`).
4. `typing_report.py` wrote `output/typing_20261030/typing_report_dg`.

Results:
- No acceptance assertion was broken.
- 3,199,606 cells, of which 509,408 have a label different from the v19 typing (option 1 on 10-30: 508,972).
- Endothelium Lymphatic: 2,705 cells (option 1: 3,033). Every other label differs from option 1 by fewer than
  150 cells.

**v20 delivery on HPG** (owner's go-ahead; job 48777105; clive, 24 cores, 300 GB; 28 min):
- **Command.** `run_v20.py --phase deliver` (`runs/20261101T165721.675277Z`), with
  `RUN_V20_FULL_QUPATH_WORKERS=8`. The variable, added for HPG allocations, defaults to 20; each image is imported
  independently.
- **Timings.** Assembly 691 s, giving `data/development/uniform-20donor-20261030-v20`, with `git_commit` 9b1da13
  recorded in `uniform_release.json`. The full QuPath project took 268 s
  (`data/development/qupath-review-20261030-v20`, `output/typing_20261030/full_qupath_assignments_v20`). The final
  reconciliation took 435 s (`output/typing_20261030/final_reconciliation_v20`).
- **Reconciliation.** Status `pass` with no errors, 20 of 20 images `pass`. Of 3,525,977 raw detections, 3,199,606
  are QC-eligible, typed, assembled and in the review project; 326,371 are excluded.
- **`raw_review_status`.** This is carried from the v19 QC reviews of the raw augmentation geometry
  (`final_reconciliation.py:407-412`):
  - 16 images: `not_reviewed_superseded_before_consumption`;
  - Aab_6450 and T1D_6533: `needs_review`;
  - ND_0112: `quality_hold`;
  - Aab_6505: `acceptable_with_notes`.

**Return to the workstation** (`output/hpg_transfer_20261031/return.sh`):
- Brought back, after confirming that none of these existed on the workstation:
  - the v20 root;
  - the v20 QuPath project;
  - `assembly_v20`, `full_qupath_assignments_v20`, `final_reconciliation_v20`;
  - `typing_report_dg`, `lymphatic_despeckle_gate_v1`;
  - the five 11-01 run journals;
  - `output/hpg_equivalence_20261031`;
  - the seal and the HPG gate dump.
- `output/hpg_return_20261101/SHA256SUMS`: 5,532 files, made on HPG (job 48780847) and verified on the workstation.
- All 20 image URIs of the v20 QuPath project and all 20 image links of the v20 root resolve on the workstation.
- HPG remains the site of record for the working run `uniform-20donor-20261029-v19-qc-work`: its DG binding, its
  structures and typing outputs, and its journals. The workstation copy is the state captured before the move.
