# Withdrawn SAM R5 boundary candidate

These files record the September 17 endocrine-seeded SAM candidate. R5 was
withdrawn after repeated HPAP-078 background-overreach errors. Its masks and
green/yellow review labels are historical evidence, not accepted islet
boundaries. See [the method and withdrawal record](../../docs/islet_regeneration_20260917.md).

The inference, finalization, gallery rendering, and QuPath import files were
archived from `scripts/cohort_expansion/`. The mask-geometry test was archived
from `scripts/tests/`. The historical preparation helper remains in
`scripts/cohort_expansion/` because later boundary diagnostics import its
shared functions. The package helper also remains there and is still covered
by the withdrawal test. The prompt-placement and withdrawal tests remain in
`scripts/tests/`.

The reproduction commands in the method record use these archive paths. They
must not be treated as the current six-donor boundary workflow or as an app
promotion recipe. The archived Python scripts import the shared preparation
and image helpers from the repository's `scripts/cohort_expansion/` directory
only when run explicitly.
