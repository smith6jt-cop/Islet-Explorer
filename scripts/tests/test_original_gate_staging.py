"""Small checks for preserving missingness and marker order in gate staging."""

from pathlib import Path
import sys

import h5py
import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from stage_original_gate_calls import frozen_gates, stage_donor  # noqa: E402


def test_original_gates_follow_marker_names_and_leave_missing_calls_null(tmp_path):
    markers = [f"m{i}" for i in range(31)]
    gates_path = tmp_path / "frozen.h5ad"
    with h5py.File(gates_path, "w") as h5:
        gate_group = h5.create_group("uns/gates")
        gate_group.create_dataset("markers", data=np.asarray(markers, dtype="S"))
        for donor in range(15):
            gate_group.create_dataset(str(donor), data=np.ones(31))
    recovered_markers, gates = frozen_gates(gates_path)
    assert recovered_markers == markers

    source = tmp_path / "78.h5ad"
    reversed_markers = list(reversed(markers))
    raw = np.zeros((3, 31), dtype=np.float32)
    raw[0, reversed_markers.index("m0")] = np.expm1(1.5)
    raw[1, :] = np.nan
    raw[2, reversed_markers.index("m1")] = np.expm1(2.5)
    with h5py.File(source, "w") as h5:
        h5.create_dataset("var/marker", data=np.asarray(reversed_markers, dtype="S"))
        h5.create_dataset("obs/cell_id", data=np.asarray(["a", "b", "c"], dtype="S"))
        h5.create_dataset("X", data=raw)
        h5.create_dataset("layers/log", data=np.log1p(raw))
    output = tmp_path / "calls.parquet"
    result = stage_donor(source, output, recovered_markers, gates, chunk_size=2)
    calls = pq.read_table(output, columns=["cell_id", "m0", "m1"]).to_pydict()
    assert calls == {"cell_id": ["a", "b", "c"],
                     "m0": [True, None, False], "m1": [False, None, True]}
    assert result["measured_cells"] == 2
    assert result["wholly_unmeasured_cells"] == 1

