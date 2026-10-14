import sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_protein_embedding import protein_pca, fit_embedding


def test_scaling_rank_and_constants():
    x = np.array([[1, 2, 9], [2, 4, 9], [4, 8, 9], [8, 16, 9.]])
    scores, params = protein_pca(x)
    assert scores.shape == (4, 1)
    assert params['active'].tolist() == [True, True, False]
    assert np.all(params['loadings'][2] == 0)
    for bad in (np.ones((5, 3)), np.array([[1., np.nan]] * 5)):
        with pytest.raises(ValueError):
            protein_pca(bad)


def test_embedding_preserves_measurements_and_ignores_scvi():
    rng = np.random.default_rng(12)
    source = ad.AnnData(rng.uniform(size=(48, 6)),
        obs=pd.DataFrame({'imageid': np.repeat(['1','2','3'], 16),
                          'count_Beta cell': np.arange(48)}, index=[f'i{x}' for x in range(48)]))
    source.layers['log'] = source.X.copy()
    result, report = fit_embedding(source, stability=False)
    np.testing.assert_array_equal(source.X, result.X)
    pd.testing.assert_series_equal(source.obs['count_Beta cell'], result.obs['count_Beta cell'])
    assert result.uns['neighbors']['params']['use_rep'] == 'X_protein_pca'
    assert report['scvi_required'] is False
    source.obsm['X_scVI_mean'] = rng.normal(size=(48, 10)) * 1e6
    second, _ = fit_embedding(source, stability=False)
    np.testing.assert_array_equal(result.obsm['X_umap'], second.obsm['X_umap'])
    assert result.obs['leiden_0.5'].equals(second.obs['leiden_0.5'])
    assert 'leiden_0.5' not in source.obs
