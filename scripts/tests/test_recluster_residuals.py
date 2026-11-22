import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "uniform"))
from recluster_residuals import constraints


def test_cluster_similarity_requires_phenotype_evidence():
    markers = ["CD45", "CD8a", "CD4", "CD3e", "CK19", "ECAD", "panCK", "INS", "CD31"]
    panel = pd.DataFrame({f"pos__{m}": [0] * 4 for m in markers})
    panel["start_class"] = ["Tissue", "Tissue", "EndocrineCluster", "Tissue"]
    panel["n_positive"] = [1, 2, 3, 2]
    panel["fold__CD4"] = 0.0
    panel["fold__CD8a"] = 2.0
    panel["in_Blood vessel"] = [False, False, False, True]
    panel["crossed_frac_vessel"] = [0.0, 0.0, 0.0, 0.8]
    panel["pos__CD45"] = 1
    panel.loc[[1, 2], "pos__CD8a"] = 1
    panel.loc[2, "pos__INS"] = 1
    panel.loc[[2, 3], "pos__CD31"] = 1
    labels = np.array(
        ["Tissue", "CD8a Tcell", "Acinar cell", "Beta cell", "Endothelium Blood Vessel"]
    )
    mask = constraints(panel, labels)
    assert list(labels[mask[0]]) == ["Tissue"]
    assert list(labels[mask[1]]) == ["Tissue", "CD8a Tcell"]
    assert list(labels[mask[2]]) == ["Tissue", "CD8a Tcell", "Beta cell"]
    assert list(labels[mask[3]]) == ["Tissue", "Endothelium Blood Vessel"]
