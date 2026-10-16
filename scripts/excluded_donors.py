"""Donors removed from the app and all analysis by owner direction.

`excluded_donors.csv` is the single record. Cohort loaders call
`assert_no_excluded()` so an excluded donor stops a run instead of being
filtered silently. Identifiers are compared after dropping an HPAP/HDL prefix,
separators and leading zeros, so 78, 0078, HPAP-078 and HDL078 all match.
The Shiny app keeps the same IDs in `ISLET_EXCLUDED_DONORS` (00_globals.R).
"""
from __future__ import annotations

import csv
from pathlib import Path
import re

REGISTRY = Path(__file__).resolve().with_name('excluded_donors.csv')
_ID = re.compile(r'(?i)^\s*(?:HPAP|HDL)?[-_ ]?0*(\d+)(?:\.0+)?\s*$')


def canonical_id(value) -> str:
    """Return the bare numeric donor ID used for comparison."""
    match = _ID.match(str(value))
    return match.group(1) if match else str(value).strip()


def excluded_ids(path=REGISTRY) -> set[str]:
    """Canonical IDs of every excluded donor, including listed aliases."""
    with Path(path).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    ids = set()
    for row in rows:
        ids.add(canonical_id(row['case_id']))
        ids.update(canonical_id(alias) for alias in row['aliases'].split('|') if alias.strip())
    return ids


def find_excluded(values, path=REGISTRY) -> list[str]:
    """The given identifiers that name an excluded donor, in input order."""
    blocked = excluded_ids(path)
    return [str(v) for v in values if canonical_id(v) in blocked]


def assert_no_excluded(values, context, path=REGISTRY):
    """Raise if any identifier names an excluded donor."""
    found = find_excluded(values, path)
    if found:
        raise ValueError(f'{context}: excluded donor(s) {sorted(set(found))} present; '
                         f'see {Path(path).name} (removed from the app and all analysis)')
