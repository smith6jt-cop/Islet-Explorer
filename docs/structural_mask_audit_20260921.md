# Structural image audit, 2026-09-21

The schema-4 candidate contains explicit structural-mask dispositions for all
21 donors. Availability differs by measurement; a donor remains discoverable
in the app even when one structural endpoint is unavailable. These are current
source-image audits, not reused acceptance from the original release.

The authoritative generated records are
`data/development/all-donor-20260920/structural/image_audit_20260921.json`
and `structural/{case_id}.json`. The audit registry binds all 126 decisions
to the exact mask and review-image hashes. It records operational sampled
image acceptance, not an independent pixel-level accuracy estimate.

## Methods and evidence

`scripts/research/structural_masks.py` reads the corrected candidate OME channel
names. Analysis geometry uses the source QuPath calibration of 0.3774 µm per
pixel; the acquired OME calibration is recorded separately. Confirmed technical
exclusions are omitted from reference sampling. Downstream spatial measurement
must intersect masks with the final assessable compartment geometry.

New ColIV and duct masks were proposed from sampled nonzero source pixels using
two-component log1p intensity mixtures, with zero intensity treated separately.
Poorly separated fits used a conservative quantile proposal. None of those
statistical fits automatically validated a mask. No disease label, endocrine
outcome or hypothesis result selected a threshold.

ColIV proposals preserve the actual threshold-positive pixels. Duct-positive
pixels require CK19 support and either panCK support or stronger CK19 signal.
A separately recorded duct envelope closes small gaps and fills lumina for
distance measurement. The envelope is never substituted for positive-pixel
area. Pixel polygonization preserves holes and exact pixel area. QuPath
annotations contain names and classification colors.

The inherited Nerve, Capillary and Lymphatic outlines were retained as proposals
and reviewed against the current correctly named channels. Every donor received
four source-tissue review fields per structure, selected without disease labels
across sampled intensity quantiles. Companion panels show:

| Structure | Source channels examined |
| --- | --- |
| ColIV | ColIV, ECAD, panCK |
| Duct | CK19, panCK, ECAD, BCatenin |
| Nerve | B3TUBB, PGP9.5, GAP43 |
| Capillary | CD31, CD34, SMA |
| Lymphatic | PDPN, CD31, HLADR |

The three inherited structures also received median, 90th-percentile and largest
assessable object-area review fields. These object-targeted fields complement
the tissue samples, where rare structures can be absent. All fields were
actually inspected. Reference locations, source identities and image hashes
remain in each donor's proposal record. Additional exact-cell source-image
reviews can withdraw an earlier acceptance, as happened for Nerve/6533.

## Current availability

| Measurement | Accepted donors | Held donors and reason |
| --- | ---: | --- |
| ColIV positive pixels | 0/21 | All donors: proposed threshold masks either captured diffuse speckle or missed substantial visible matrix. No ColIV classifier or H7 endpoint is released. |
| Duct positive pixels | 18/21 | 6533, 6538, 6550: substantial CK19/panCK-positive structures are missed. |
| Duct envelope for distance | 18/21 | Same three donors; an incomplete positive mask cannot certify its derived envelope. |
| Nerve area | 18/21 | 6505, 6520, 6533: inherited outlines include conspicuous technical signal. |
| Capillary area | 11/21 | 112, 6450, 6479, 6505, 6516, 6521, 6533, 6534, 6548, 6551: selected masks include walls of large SMA-coated vessels and are not certified as capillary-specific area. |
| Lymphatic area | 1/21 | Only 6520 has accepted sampled continuous PDPN vessel profiles; the remaining masks mix plausible vessels with puncta, diffuse signal or nonspecific bright objects. |

The large SMA-ring vessel in donor 6090 was inspected at full resolution and is
excluded from its selected Capillary mask; adjacent narrow CD31 profiles are
selected. Its capillary acceptance is retained.

Technical findings requiring explicit regional exclusion or mask reassessment:

- 6505: saturated ring and broad halo, near source crop `(776, 1857, 256, 256)`,
  inherited Nerve feature 7439; similar signal enters the Lymphatic mask.
- 6520: diffuse saturated edge, near `(11083, 5556, 256, 256)`, inherited Nerve
  feature 5962. An additional tissue review field starts at `(11162, 5286)`.
- 6533: B3TUBB source review centered at `(5213.418, 862.091)` reveals a broad
  saturated area with rectangular intensity boundaries. The inherited Nerve
  mask covers 12,743.38 of 16,384 pixels in the surrounding 128×128-pixel box.
  The additional source panel and UUID evidence are bound in the registry;
  its initial sampled acceptance was withdrawn.

Held masks remain available as review artifacts. They yield missing structural
measurements, not zero density. Empty masks are not declared assessably empty
without explicit evidence. Biological PDPN or B3TUBB positivity, cell identity
and structural-mask specificity are separate decisions.

## Independent marker-positivity review

`positivity_restore_v1/nonimmune_review_20260921.json` contains 252 explicit
donor-marker decisions for CD31, CD34, ECAD, panCK, CK19, KRT14, SMA, VIM,
PGP9.5, B3TUBB, GAP43 and PDPN. All 36 source-image pages were reviewed, each
showing seven donor strips with reference-control, threshold-adjacent and
high-target examples. Five additional epithelial donor strips were reopened.
The calls are bound to current fresh RESTORE evaluations, native cell UUIDs,
source panels and their hashes. A computational failure cannot be overruled
by visual review. This partial registry must be merged with all remaining
marker dispositions before producing the complete 31-marker call tables.

Accepted donor counts are CD31 21, CD34 21, ECAD 20, B3TUBB 16, CK19 15,
SMA 14, VIM 14, PGP9.5 12, panCK 11, GAP43 10, PDPN 5 and KRT14 2.
Dim or poorly attributed examples remain unresolved; this is not a negative
biological call. These sampled operational decisions do not establish
independent sensitivity or specificity, or assign a lineage by one marker.

## Verification

Focused tests check exact pixel area and holes, the distinction between duct
positive area and envelope, degenerate-marker handling, missing companion
evidence, modified masks and modified source review images. The acceptance
tool validates all decisions before writing any donor registry. Current
structural and schema-4 Python checks pass (12 tests). App-specific schema-4
tests and isolated browser evidence are separate from scientific acceptance.
