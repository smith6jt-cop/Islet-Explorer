# Deploying Islet Explorer to UF Research Computing PubApps

Runbook for hosting the app **public-facing** on HiPerGator's **PubApps** platform
(open internet + RC-managed SSL), replacing the on-prem lab deployment
(`http://10.15.152.7:8080/islet-explorer/`).

- **Access level:** Public Open (no login).
- **Runtime:** Apptainer/Podman container (RC does **not** allow Docker).
- **Includes** the ~109 GB OME-TIFF Viewer image store.

Full rationale and the app inventory that drives these choices are in the plan file:
`~/.claude/plans/review-the-repo-code-reflective-hanrahan.md`.

## Files in this directory

| File | Purpose |
|---|---|
| `islet-explorer.def` | Apptainer image recipe (OS libs + R packages + anndata venv + shiny-server). |
| `shiny-server.conf` | Serves the app at `/islet-explorer/` on `:3838`; baked into the image. |
| `run.sh` | Bind mounts + secret injection to launch the container on the VM. |
| `Renviron.pubapps.example` | Template for the deployment `.Renviron` (no secrets). |

`scripts/install_shiny_deps.R` (repo root) is the authoritative R package manifest the image builds from.

---

## ⚠️ Before anything leaves the lab

The current `app/shiny_app/.Renviron` holds a plaintext UF Navigator API key
(`KEY=sk-...`). **Rotate/revoke it** and never place the real key on the shared `/pubapps`
volume or in the image. It is injected at runtime via `run.sh`'s `--env-file` instead.

## Step 0 — Provision (UFRC)

1. Buy a **PubApps** instance ($300/yr) via the HiPerGator Service Purchase form (or a Trial
   Allocation to pilot). Base = 2 vCPU / 16 GB RAM / **1 TB `/pubapps`**.
2. Open a ticket at https://support.rc.ufl.edu with the checklist:
   - Access: **Public Open**. Runtime: **Apptainer/Podman**. App: R **shiny-server**, port **3838**, path **`/islet-explorer/`**.
   - Domain: **`islet-explorer.rc.ufl.edu`** (RC SSL).
   - Reverse proxy **must** support **WebSocket upgrade** (Shiny) **and unbuffered HTTP Range** (Viewer streams multi-GB OME-TIFFs).
   - **Confirm `/pubapps` allocation is 1 TB** (the oft-quoted 40 GB is the *home* quota, a different filesystem — not used here). ~111 GB of data+images fits 1 TB with margin.
   - No SLURM jobs; persistent process. Needs outbound HTTPS to `api.ai.it.ufl.edu`.
   - Give developer GatorLink usernames for SSH setup.
3. RC returns: VM hostname, `$serviceuser`, `/pubapps/$project`, `/podman/$serviceuser`, public URL mapping.

## Step 1 — Build the image (on a HiPerGator node with Apptainer)

```bash
cd <repo root>
# --fakeroot lets an unprivileged HiPerGator user run the %post apt-get/install steps.
apptainer build --fakeroot deploy/islet-explorer.sif deploy/islet-explorer.def   # %files paths need CWD=repo root
# copy to the VM's container store:
rsync -av deploy/islet-explorer.sif $VM:/podman/$SERVICEUSER/
```

The image installs only the **required** R manifest. `scripts/install_shiny_deps.R` now exits
non-zero if any required package fails, so a broken build fails here rather than as an HTTP 500
from a shiny-server worker. The optional stacks are opt-in and none are needed by the app (the
trajectory ships precomputed in the H5AD):

```bash
# only if you need the offline trajectory/recompute paths:
INSTALL_BIOC=1 INSTALL_VITESSCER=1 INSTALL_LAMIAN=1 Rscript scripts/install_shiny_deps.R
```

## Step 2 — Stage code + data + images on `/pubapps`

Target layout (the app's `../../data` requires `data/` to be the app dir's grandparent):

```
/pubapps/$project/islet-explorer/
  app/shiny_app/                 # git checkout (code is tracked)
    www/local_images/            # 109 GB OME-TIFFs (gitignored — rsync separately)
    .Renviron                    # from deploy/Renviron.pubapps.example (NO secret)
  data/                          # 1.9 GB (gitignored — rsync separately)
```

```bash
# on the VM (or via HiPerGator, since /pubapps isn't visible from HPG compute):
git clone <repo> /pubapps/$project/islet-explorer

# stage the gitignored trees from the lab machine (long — use tmux + --partial):
rsync -av --progress --partial data/ \
      $VM:/pubapps/$project/islet-explorer/data/
rsync -av --progress --partial app/shiny_app/www/local_images/ \
      $VM:/pubapps/$project/islet-explorer/app/shiny_app/www/local_images/

# deployment .Renviron (no secret) + the runtime secret (chmod 600, NOT in git):
cp deploy/Renviron.pubapps.example \
   /pubapps/$project/islet-explorer/app/shiny_app/.Renviron
mkdir -p /pubapps/$project/secrets
printf 'KEY=%s\n' 'sk-...ROTATED...' > /pubapps/$project/secrets/islet.env
chmod 600 /pubapps/$project/secrets/islet.env
```

## Step 3 — Launch

```bash
PROJECT=$project SERVICEUSER=$serviceuser bash deploy/run.sh
```
RC then points the reverse proxy at `127.0.0.1:3838` and maps the public host/SSL.
For production, RC wires this as a persistent systemd/quadlet service (see `run.sh` for the
exact mounts/env it must replicate). **Updates:** `git pull` in `/pubapps/.../islet-explorer`
+ `touch app/shiny_app/restart.txt` — no rebuild unless dependencies change.

## Step 4 — Verify (end-to-end)

Open `https://islet-explorer.rc.ufl.edu/islet-explorer/` and check:

1. App renders; WebSocket stays connected (no "disconnected from server").
2. **Plot → Population Measure lists all ~21 cell types** — the key test that H5AD loaded via
   reticulate/anndata rather than the Excel fallback.
3. **Viewer** renders an OME-TIFF (Range serving + correct public URL through the proxy).
4. **AI Assistant** returns a reply (KEY injected + egress to `api.ai.it.ufl.edu`).
5. **Spatial Card C** distance views load (`cell_distance.csv` staged/readable).
6. Container/shiny-server logs show no anndata/reticulate errors. `run.sh` bind-mounts
   `/pubapps/$PROJECT/logs` → `/var/log/shiny-server`, so read them on the host with
   `tail -f /pubapps/$PROJECT/logs/*.log`. Without that bind `--writable-tmpfs` discards
   them — which is precisely when you need them (a worker that dies at startup leaves nothing).

## Resource sizing

Base 2 vCPU / 16 GB / 1 TB is adequate for a public research showcase at light/moderate
concurrency: ~1.5–2.5 GB per Shiny worker (→ ~6–8 workers in 16 GB); 2 cores is the limiting
factor under many simultaneous users (Shiny is single-threaded per session). Scale later via
**NCU** (CPU/RAM) and **BlSU** (storage) purchases + a VM restart — no re-architecture.

## Open items to confirm with RC

- Exact public path (`/` vs `/islet-explorer/`) and whether they rewrite it (the app handles
  either via `session$clientData` in `R/viewer_helpers.R`, but the mapping must be agreed).
- WebSocket + unbuffered Range support on the proxy.
- Rootless `run_as`/UID mapping so the shiny worker can read `$serviceuser`-owned `/pubapps` files
  (see the comment in `shiny-server.conf` and `--userns=keep-id` note in `run.sh`).
