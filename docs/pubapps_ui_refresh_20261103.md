# PubApps UI refresh build — 2026-11-03

## Final Production Release — 2026-11-04

The earlier local-review candidate described below was superseded after the
parity verifier was updated for the Unclassified selector helpers. The final
release `pubapps-2026-11-ui-refresh-1659c9d` is active on production. The
verified 20-donor schema-5 bundle passed **25/25** parity checks against the
canonical component-app-v2 source. The validator commit `7ffa7a1` is recorded
separately from producer/image commit `1659c9d`; the verifier accepts the
validator-only commit difference only when its sole changed file is the
parity verifier.

| Item | Final identity or result |
|---|---|
| Application/data producer | `1659c9dd01d8d409030cfb4618eb9f3edf56488e` |
| Parity validator | `7ffa7a1f38a95a23f311e62892df857f2f16b466` |
| OCI image | `sha256:16c4979fff4d7e8d72befbaade013ce668153dc6161f938c8bc6a49a60e5ae24` |
| Release manifest | `a097742540f8bbcc7bb532ecd031dc7f27b650dd03613a3148ef0abca710d50b` |
| Data manifest | `fd21101e781d1b713d3f5a26f65cf2249e7591240498cd5238510bd2608e2f15` |
| Release parity | 25 passed, 0 failed; exact categorical/count/NA parity and numeric tolerance `1e-12` |
| Host verification | Full hashes, schema 5, image identity, read-only scientific mounts, and TIFF Range response passed |
| Public proxy | Authenticated reviewer access passed 15/15 checks |
| Activated | `2026-11-04T13:18:01Z` on `pubmaigan1`, production port 3838 |
| Previous release | `pubapps-2026-11-components-v1`, retained for rollback |

The production data has 20 TIFFs plus 4,907 other data files. The image tree
is byte-identical to the prior release and is hard-linked on the PubApps
filesystem; no TIFF bytes were copied to the host. A sampled TIFF has the same
device and inode in both release directories. The app and read-only image mount
use their normal paths. The dev service was not changed, and the local preview
has been stopped. The clean off-campus acceptance run, reviewer credential
delivery, and VM reboot persistence test remain pending. See the
[deployment receipt](evidence/pubapps_live_20261104_ui_refresh.json),
[25-check parity report](evidence/release_parity_pubapps_ui_refresh_20261104.json),
and [reviewer-access report](evidence/reviewer_access_pubapps_ui_refresh_20261104.json).

The original `26b6613` local-review build and no-deploy statement below are
historical and describe the state before the final release was rebuilt and
activated.

The PubApps release candidate `pubapps-2026-11-ui-refresh-26b6613` is built
from the exact source and image running in the local preview. It pairs those
UI changes with the unchanged schema-5 data bundle from the accepted production
release `pubapps-2026-11-components-v1`. No new scientific data was introduced.

| Item | Identity |
|---|---|
| Application source | `26b66131668c01326d4b95bdb1efbf1d37ef7527` |
| OCI image | `sha256:2931955ab95330c7804bfea7f9ca25774aa4813f7c9db645a3fd63e35f99035e` |
| Image security | Strict gate passed; 0 High, 0 Critical, no DuckDB deferral |
| Data producer | `67ee6c9d0b28f467eee589fea24829a6546e3d74` |
| Data manifest | `8ddeecf3639921d2e53538e2d86b6bb4278cc622d29a928a3136fdece3d9f6b1` |
| Release manifest | `8c43f6823bd52dddc51efa0912541751122987326cda5453cc2e0b0b1faf9d6e` |
| Archive | `tmp/pubapps-ui-refresh-20261103/islet-explorer-pubapps-2026-11-ui-refresh-26b6613.oci.tar` |
| Archive SHA-256 | `a3a73bae10534f1507c37574647d0169a7ff7886e74a292d881eda6cac215874` |

The release manifest and full file verification passed for 20 OME-TIFFs and
4,907 non-image data files. The image is the strict-gate-passed artifact already
built for the local preview; its source revision matches the packaged release.
The OCI archive and checksum metadata are ready for transfer. This is a local
build only: it has not been transferred to HiPerGator, staged on PubApps, or
activated on dev or production. The live services were not changed.

Machine-readable identities and checks are recorded in
[`docs/evidence/pubapps_ui_refresh_20261103.json`](evidence/pubapps_ui_refresh_20261103.json).
