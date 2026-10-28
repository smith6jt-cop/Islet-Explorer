"""registry.py -- marker roles and typing constants for one panel and one tissue, read from data, never from code.

Owner, 2026-10-28: "The lineage and shared lists are not static and will be panel and tissue specific which needs to be
clear for the pipeline as a whole." Every typing module asks this module which markers it may use and how; no marker
list is written into a typing module. A run names its panel and tissue in `scripts/uniform/cohort.csv` (columns
`panel`, `tissue`), and the files are

    registry/<panel>__<tissue>.csv              one row per channel: role, compartment, lineage, measurement, ...
    registry/<panel>__<tissue>__constants.csv   one row per constant: value, source, reason

Roles (canonical marker biology for the tissue, owner 2026-10-28):
  defining     a positive value places a cell in the marker's compartment (and, inside it, the lineage)
  shared       canonically on several compartments; adds evidence inside a compartment a defining marker found,
               never places a cell
  state        measured and flagged (e.g. Ki67), never part of an identity
  none         acquired but not used in typing (e.g. ColIV, owner 2026-10-28: "excluded from typing completely")
  never_used   not acquired on every image of the run (owner 2026-10-28: "Don't use markers that are not in the panel
               for every image")
Measurement: one REDSEA compartment per marker for every cell, whole cell, membrane or nucleus (owner, 2026-10-28:
"We can't have two ways of typing. We need to use whole cell and/or membrane and/or nucleus for all measurements.").
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
COHORT = HERE.parent / 'cohort.csv'
REGISTRY_DIR = HERE / 'registry'
TYPING_ROLES = ('defining', 'shared')
COMPARTMENT_COLUMN = {'whole_cell': 'Cell', 'membrane': 'Membrane', 'nucleus': 'Nucleus'}


def split(s: str) -> list[str]:
    return [x for x in (s or '').split(';') if x]


@dataclass
class Registry:
    name: str
    table: pd.DataFrame
    constants: dict = field(default_factory=dict)
    constant_sources: dict = field(default_factory=dict)

    # -- markers ----------------------------------------------------------------------------------------------------
    def role(self, marker: str) -> str:
        return self.table.loc[marker, 'role']

    def markers(self, *roles: str) -> list[str]:
        roles = roles or TYPING_ROLES
        return [m for m in self.table.index if self.table.loc[m, 'role'] in roles]

    def typing_markers(self) -> list[str]:
        """Markers with a typing role (defining or shared)."""
        return self.markers(*TYPING_ROLES)

    def measured_markers(self) -> list[str]:
        """Markers whose positivity is computed: typing markers and states."""
        return self.markers('defining', 'shared', 'state')

    def compartments(self) -> list[str]:
        return list(dict.fromkeys(self.table.loc[self.markers('defining'), 'compartment']))

    def defining(self, compartment: str) -> list[str]:
        t = self.table
        return [m for m in self.markers('defining') if t.loc[m, 'compartment'] == compartment]

    def shared(self, compartment: str) -> list[str]:
        t = self.table
        return [m for m in self.markers('shared') if compartment in split(t.loc[m, 'shared_by'])]

    def lineages(self, marker: str) -> list[str]:
        return split(self.table.loc[marker, 'lineage'])

    def positivity_class(self, marker: str) -> str:
        return self.table.loc[marker, 'positivity_class']

    def measurement(self, marker: str) -> str:
        return self.table.loc[marker, 'measurement']

    def candidates(self, marker: str) -> list[str]:
        return split(self.table.loc[marker, 'measurement_candidates']) or [self.measurement(marker)]

    def column(self, marker: str, measurement: str | None = None) -> str:
        """The REDSEA column holding this marker's measurement, e.g. 'Cell__INS'."""
        return f"{COMPARTMENT_COLUMN[measurement or self.measurement(marker)]}__{marker}"

    # -- constants --------------------------------------------------------------------------------------------------
    def const(self, name: str) -> float:
        if name not in self.constants:
            raise KeyError(f'{self.name}: constant {name} is not in the constants file')
        return self.constants[name]

    # -- checks -----------------------------------------------------------------------------------------------------
    def validate(self, image_columns: dict[str, set[str]]) -> None:
        """Fail if a marker with a role is missing on any image, or a never_used marker is marked present everywhere."""
        problems = []
        for m in self.measured_markers():
            for img, cols in image_columns.items():
                if self.column(m) not in cols:
                    problems.append(f'{m} ({self.column(m)}) missing on {img}')
            if self.table.loc[m, 'on_every_image'] != 'yes':
                problems.append(f'{m} has a role but is not acquired on every image')
        if problems:
            raise SystemExit(f'{self.name}: ' + '; '.join(problems))


def run_panel_tissue(cohort: Path = COHORT) -> tuple[str, str]:
    rows = list(csv.DictReader(cohort.open()))
    pairs = {(r.get('panel', ''), r.get('tissue', '')) for r in rows}
    if len(pairs) != 1 or '' in next(iter(pairs)):
        raise SystemExit(f'{cohort}: every row needs the same non-empty panel and tissue, found {sorted(pairs)}')
    return next(iter(pairs))


def load(panel: str | None = None, tissue: str | None = None, directory: Path = REGISTRY_DIR) -> Registry:
    if panel is None or tissue is None:
        panel, tissue = run_panel_tissue()
    name = f'{panel}__{tissue}'
    table = pd.read_csv(directory / f'{name}.csv', keep_default_na=False, dtype=str).set_index('marker')
    consts, sources = {}, {}
    cpath = directory / f'{name}__constants.csv'
    if cpath.exists():
        for r in csv.DictReader(cpath.open()):
            consts[r['name']] = float(r['value']) if r['value'] not in ('', 'pending') else float('nan')
            sources[r['name']] = r['source']
    return Registry(name, table, consts, sources)
