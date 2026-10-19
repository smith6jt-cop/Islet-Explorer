import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
from composition_embedding import ENDOCRINE, cross_donor_disease, endocrine_shares, fit_composition_embedding

NAMES = list(ENDOCRINE) + ['Exocrine']


def source(n=160, seed=3):
    rng = np.random.default_rng(seed)
    counts = rng.integers(0, 12, size=(n, len(NAMES))).astype(np.float32)
    counts[:, NAMES.index('Endocrine')] = 0      # every endocrine cell has a subtype, as on v5
    counts[0, :len(ENDOCRINE)] = 0                # an object with no endocrine cell in its core
    core = counts.sum(axis=1)
    donors = np.repeat(['1', '2', '3', '4'], n // 4)
    status = np.repeat(['ND', 'Aab+', 'T1D', 'ND'], n // 4)
    a = ad.AnnData(rng.uniform(size=(n, 3)).astype(np.float32),
                   obs=pd.DataFrame({'imageid': donors, 'donor_status': status, 'total_cells_core': core,
                                     'object_class': 'Islet',
                                     'leiden_0.5': pd.Categorical(rng.integers(0, 3, n).astype(str))},
                                    index=[f'o{i}' for i in range(n)]))
    a.obsm['phenotype_counts'] = counts
    a.obsm['X_umap'] = rng.normal(size=(n, 2))
    a.uns['phenotype_names'] = np.asarray(NAMES, dtype=object)
    a.uns['protein_embedding'] = {'donor_correction': 'harmony'}
    return a


def test_shares_are_of_core_endocrine_cells():
    a = source()
    shares = endocrine_shares(a)
    assert np.isnan(shares[0]).all()
    np.testing.assert_allclose(shares[1:].sum(axis=1), 1)
    assert (shares[1:, ENDOCRINE.index('Endocrine')] == 0).all()


def test_fit_covers_only_ten_cell_objects_and_keeps_the_protein_record():
    a = source()
    before = a.copy()
    result, report = fit_composition_embedding(a, stability=False)
    shares = endocrine_shares(a)
    fit = (a.obs.total_cells_core.to_numpy() >= 10) & np.isfinite(shares).all(axis=1)
    assert report['objects'] == fit.sum() and report['objects_left_out'] == (~fit).sum()
    assert np.isfinite(result.obsm['X_umap'][fit]).all() and np.isnan(result.obsm['X_umap'][~fit]).all()
    np.testing.assert_array_equal(result.obs.leiden_umap_1.to_numpy(), result.obsm['X_umap'][:, 0])
    assert result.obs['leiden_0.5'][fit].notna().all() and result.obs['leiden_0.5'][~fit].isna().all()
    assert result.obs.in_composition_embedding.tolist() == fit.tolist()
    # protein results move under their own names; the source is untouched
    np.testing.assert_array_equal(result.obsm['X_umap_protein'], before.obsm['X_umap'])
    assert result.obs['protein_leiden_0.5'].equals(before.obs['leiden_0.5'])
    assert 'role' in result.uns['protein_embedding'] and 'role' not in a.uns['protein_embedding']
    np.testing.assert_array_equal(a.obsm['X_umap'], before.obsm['X_umap'])
    method = result.uns['composition_embedding']
    assert method['donor_correction'] == 'none' and method['min_core_cells'] == 10
    assert list(method['excluded_constant_features']) == ['Endocrine']
    again, _ = fit_composition_embedding(a, stability=False)
    np.testing.assert_array_equal(result.obsm['X_umap'], again.obsm['X_umap'])


def test_too_few_objects_is_an_error():
    a = source()
    a.obs['total_cells_core'] = 5.0
    with pytest.raises(ValueError):
        fit_composition_embedding(a, stability=False)


def test_cross_donor_disease_ignores_same_donor_neighbours():
    # Two donors per status, each donor a tight blob, the two same-status donors far apart: every neighbour from
    # another donor has the other disease.
    rng = np.random.default_rng(0)
    centres = {'a': (0, 0), 'b': (0, 1), 'c': (50, 0), 'd': (50, 1)}
    donor = np.repeat(list(centres), 40)
    status = np.where(np.isin(donor, ['a', 'c']), 'ND', 'T1D')
    coords = np.vstack([np.asarray(centres[d]) + rng.normal(scale=0.05, size=2) for d in donor])
    out = cross_donor_disease(coords, status, donor, k=45)
    assert out['other_donor_same_disease'] == 0
    assert out['same_donor_neighbour_fraction'] > 0.8
