# Typing processing on HiPerGator (owner decisions of 2026-10-31)

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
