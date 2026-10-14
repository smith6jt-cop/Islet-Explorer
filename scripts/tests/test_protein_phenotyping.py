import sys
from pathlib import Path
import h5py
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from protein_phenotyping import (QC_TOTAL_MFI, classify, derive_original_qc_threshold,
                                 qc_mask, read_rules, rescale_proteins)


def test_lineage_conflicts_and_subtypes():
    rules = read_rules()
    markers = sorted({m for r in rules for k in ['all_positive','any_positive','all_negative'] for m in r[k]})
    patterns = [['INS'], ['INS','GCG'], ['CD3e','CD8a'], ['CD3e','CD8a','CD4'],
                ['CK19'], [], ['SMA','B3TUBB'], ['HLADR'], ['CD31'], ['PDPN']]
    x = np.zeros((len(patterns), len(markers)))
    for i, pattern in enumerate(patterns):
        for m in pattern: x[i, markers.index(m)] = .8
    labels, reason = classify(x, markers)
    assert labels.tolist() == ['Beta cell', 'Endocrine', 'CD8a Tcell', 'T cell',
                              'Ductal', 'Unknown', 'Unknown', 'APCs', 'Blood Vessel', 'Lymphatic']
    assert reason[6] == 'conflicting_all_evidence'
    x[0, 0] = np.nan
    assert classify(x, markers)[1][0] == 'missing_measurements'
    with pytest.raises(ValueError, match='Missing required'):
        classify(x[:, 1:], markers[1:])


def test_qc_and_rescaling_missingness():
    x = np.array([[1000.,5000], [1,2], [np.nan,np.nan]])
    assert qc_mask(x).tolist() == [True,False,False]
    with pytest.raises(ValueError):qc_mask(np.array([[1,np.nan]]))
    raw = np.expm1(np.array([[0.,2], [2,2], [4,2]]))
    scaled, logged, params = rescale_proteins(raw, [2,3])
    np.testing.assert_allclose(scaled[:,0], [0,.5,1])
    np.testing.assert_allclose(scaled[:,1], 0)
    assert params['gate_outside_observed_range'] == [False,True]
    outside, _, _ = rescale_proteins(raw[:, :1], [5])
    assert (outside < .5).all()


def test_notebook_qc_excludes_missing_and_low_total_without_imputation():
    values = np.array([[2000.0, 4000.0], [2000.0, 1000.0],
                       [np.nan, np.nan], [6000.0, 0.0]])
    assert QC_TOTAL_MFI > 3000.0 + 1000.0  # second row is below the pinned threshold
    assert qc_mask(values).tolist() == [True, False, False, False]


def test_qc_threshold_derivation_uses_first_percentile_and_checks_count(tmp_path):
    rng = np.random.default_rng(0)
    raw = rng.uniform(0, 100, size=(1000, 4)).astype(np.float32)
    raw[:5, 1:] = 0.0  # fewer than two detected markers
    path = tmp_path / "unfiltered.h5ad"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("X", data=raw)
    totals = raw.astype(np.float64).sum(axis=1)
    expected_threshold = float(np.quantile(totals, 0.01))
    kept = int((((raw > 0).sum(axis=1) >= 2) & (totals >= expected_threshold)).sum())
    threshold, n, retained = derive_original_qc_threshold(path, expected_kept=kept, chunk=128)
    assert (threshold, n, retained) == (expected_threshold, 1000, kept)
    try:
        derive_original_qc_threshold(path, expected_kept=kept + 1, chunk=128)
    except ValueError:
        pass
    else:
        raise AssertionError("count mismatch must fail")
