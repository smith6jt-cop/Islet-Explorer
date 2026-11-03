# PubApps UI refresh build — 2026-11-03

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
