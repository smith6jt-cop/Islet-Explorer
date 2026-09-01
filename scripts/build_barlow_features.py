#!/usr/bin/env python
"""
Build the Barlow-style islet composition feature matrix.

Reproduces the feature construction from Barlow et al., eLife 2024;13:RP100535
("The extra-islet pancreas supports autoimmunity in human type 1 diabetes"),
Methods -> "Islet segmentation and pseudotime" -> "PAGA analysis":

    "For each islet, the number of each cell type inside the islet and between
     the islet edge and 20 microns beyond were extracted. To adjust for
     variation due to the islet size, the cell type counts were divided by the
     number of endocrine cells inside the islet. Data were then log-transformed."

That maps onto this project exactly: our core / peri split IS "inside the islet"
and "between the islet edge and 20 um beyond" (the `_exp20um` annotation).

Endocrine denominator
---------------------
Barlow merged alpha/beta/delta into a single 'Endocrine' type for islet
definition, and divide by "the number of endocrine cells inside the islet".
Our phenotyping keeps Alpha/Beta/Delta separate AND carries a residual
'Endocrine' class (any of INS/GCG/SST+ that did not resolve to a subtype), so
the faithful denominator here is their sum. Floored at 1 so an islet with no
resolved endocrine cell cannot produce an infinity.

What is NOT done here
---------------------
The min-cells filter is deliberately left to the consumer. Barlow filter at
>= 10 cells BEFORE computing the embedding, and the app recomputes pseudotime
whenever that threshold changes -- so the filter cannot be baked into this file
without freezing the threshold. Normalisation is per-islet and therefore
independent of which islets are kept, so it is safe to precompute here.

Output: data/barlow_islet_features.csv
    islet_id, case_id, donor_status, total_cells_core, total_cells_peri,
    n_endocrine_core,
    core_<type> x21, peri_<type> x21      log1p(count / endocrine) -- PAGA inputs
    n_core_<type> x21, n_peri_<type> x21  raw counts, so the app can build a
                                          core+peri combined ratio (log of a sum
                                          is not the sum of logs)

Usage:
    conda activate scvi-env
    python scripts/build_barlow_features.py
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(SCRIPT_DIR, '..')
CORE = os.path.join(ROOT, 'islet_analysis', 'islets_core_fixed.h5ad')
PERI = os.path.join(ROOT, 'islet_analysis', 'islets_peri_fixed.h5ad')
OUT = os.path.join(ROOT, 'data', 'barlow_islet_features.csv')

# Our analogue of Barlow's merged 'Endocrine' type.
ENDOCRINE = ['Alpha cell', 'Beta cell', 'Delta cell', 'Endocrine']


def count_frame(a, prefix):
    """21 count_<type> obs columns -> DataFrame indexed by islet_id."""
    cols = [c for c in a.obs.columns if c.startswith('count_')]
    df = a.obs[cols].copy()
    df.columns = [c[len('count_'):] for c in cols]
    df.index = a.obs.index.astype(str)
    return df.add_prefix(prefix)


def main():
    core = ad.read_h5ad(CORE)
    peri = ad.read_h5ad(PERI)
    print(f"core {core.shape}  peri {peri.shape}")

    c = count_frame(core, '')
    p = count_frame(peri, '')

    # Join on islet_id -- the two files are NOT in the same row order.
    common = c.index.intersection(p.index)
    print(f"islets in both core and peri: {len(common)} "
          f"(core {len(c)}, peri {len(p)})")
    c, p = c.loc[common], p.loc[common]

    endo_present = [e for e in ENDOCRINE if e in c.columns]
    missing = set(ENDOCRINE) - set(endo_present)
    if missing:
        print(f"  NOTE: endocrine types absent from phenotyping: {sorted(missing)}")
    n_endo = c[endo_present].sum(axis=1)
    print(f"  endocrine cells/islet: median {n_endo.median():.0f}, "
          f"zero in {(n_endo == 0).sum()} islets (floored at 1)")
    denom = n_endo.clip(lower=1)

    # Barlow: counts / endocrine cells, then log. log1p rather than log so that
    # a zero count stays 0 instead of -inf -- most cell types are absent from
    # most islets, so plain log would void the majority of the matrix.
    core_f = np.log1p(c.div(denom, axis=0)).add_prefix('core_')
    peri_f = np.log1p(p.div(denom, axis=0)).add_prefix('peri_')

    meta = pd.DataFrame({
        'islet_id': common,
        'case_id': core.obs.loc[common, 'imageid'].astype(str).values,
        'donor_status': core.obs.loc[common, 'donor_status'].astype(str).values,
        'total_cells_core': core.obs.loc[common, 'total_cells'].astype(int).values,
        'total_cells_peri': peri.obs.loc[common, 'total_cells'].astype(int).values,
        'n_endocrine_core': n_endo.values,
    }, index=common)

    # RAW counts as well as the log-ratios. The log-ratio columns are the PAGA
    # embedding inputs and must stay exactly as Barlow specify (core and peri as
    # SEPARATE features). The app's y-axis, however, needs a core+peri combined
    # ratio, and log(a) + log(b) != log(a + b) -- so the combination has to be
    # done from counts, not from the log columns.
    core_n = c.add_prefix('n_core_')
    peri_n = p.add_prefix('n_peri_')

    out = pd.concat([meta, core_f, peri_f, core_n, peri_n], axis=1)
    out.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")
    print(f"  {out.shape[0]} islets x {out.shape[1]} columns "
          f"({core_f.shape[1]} core + {peri_f.shape[1]} peri log-ratio features, "
          f"plus the same counts raw as n_core_* / n_peri_*)")
    print(f"  donor status: {dict(meta['donor_status'].value_counts())}")
    for t in (1, 10, 25):
        print(f"  islets with total_cells_core > {t}: "
              f"{int((meta['total_cells_core'] > t).sum())}")


if __name__ == '__main__':
    main()
