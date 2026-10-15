#!/usr/bin/env python3
"""Apply the frozen 15-donor numeric marker gates to six raw donor H5ADs.

This stages nullable per-cell marker calls only. It does not assign phenotypes,
fit scVI, or make an app dataset. The gates are the original ``log1p(raw MFI)``
thresholds in the frozen single-cell H5AD, not later calibrated proposals.

Run with the scVI environment, which supplies h5py, NumPy and PyArrow::

    python scripts/stage_original_gate_calls.py --output tmp/original-gate-calls
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
RAW = ROOT / "data/development/donor-additions-20260923-raw-v1/additions/raw"
DONORS = ("78", "6090", "6232", "6520", "6539", "6558")


def strings(dataset):
    return [value.decode() if isinstance(value, bytes) else str(value)
            for value in dataset[:]]


def frozen_gates(path: Path):
    with h5py.File(path) as h5:
        gates = h5["uns/gates"]
        markers = strings(gates["markers"])
        if len(markers) != 31 or len(set(markers)) != 31:
            raise ValueError("Frozen marker gate table is not the 31-marker panel")
        columns = [key for key in gates if key != "markers"]
        if len(columns) != 15:
            raise ValueError("Expected 15 original donor gate columns")
        reference = np.asarray(gates[columns[0]][:], dtype=np.float64)
        if not np.all(np.isfinite(reference)):
            raise ValueError("Frozen numeric gates contain non-finite values")
        for column in columns[1:]:
            if not np.array_equal(reference, np.asarray(gates[column][:])):
                raise ValueError(f"Gate table differs for original donor {column}")
    return markers, reference


def stage_donor(path: Path, destination: Path, markers: list[str], gates: np.ndarray,
                chunk_size: int):
    counts = np.zeros(len(markers), dtype=np.int64)
    measured_count = missing_count = 0
    with h5py.File(path) as h5:
        source_markers = strings(h5["var/marker"])
        if len(source_markers) != len(markers) or set(source_markers) != set(markers):
            raise ValueError(f"{path.name}: marker panel differs from frozen original")
        order = [source_markers.index(marker) for marker in markers]
        size = h5["X"].shape[0]
        if h5["X"].shape[1] != len(markers) or h5["layers/log"].shape != h5["X"].shape:
            raise ValueError(f"{path.name}: raw/log matrix shape mismatch")
        if h5["obs/cell_id"].shape[0] != size:
            raise ValueError(f"{path.name}: cell IDs do not match matrix rows")
        schema = pa.schema([("cell_id", pa.string())] +
                           [(marker, pa.bool_()) for marker in markers])
        with pq.ParquetWriter(destination, schema=schema, compression="zstd") as writer:
            for start in range(0, size, chunk_size):
                stop = min(start + chunk_size, size)
                raw = np.asarray(h5["X"][start:stop])[:, order]
                logged = np.asarray(h5["layers/log"][start:stop])[:, order]
                finite = np.isfinite(raw)
                if np.any(raw[finite] < 0):
                    raise ValueError(f"{path.name}: negative raw MFI")
                if not np.array_equal(finite, np.isfinite(logged)):
                    raise ValueError(f"{path.name}: raw/log missingness differs")
                if not np.allclose(logged[finite], np.log1p(raw[finite]),
                                   rtol=1e-6, atol=1e-6):
                    raise ValueError(f"{path.name}: log layer is not log1p(raw MFI)")
                complete = finite.all(axis=1)
                absent = ~finite.any(axis=1)
                if not np.all(complete | absent):
                    raise ValueError(f"{path.name}: partially measured cells need review")
                measured_count += int(complete.sum())
                missing_count += int(absent.sum())
                positive = logged >= gates
                counts += np.sum(positive & finite, axis=0)
                ids = strings(h5["obs/cell_id"][start:stop])
                if len(set(ids)) != len(ids):
                    raise ValueError(f"{path.name}: duplicate cell ID within chunk")
                arrays = [pa.array(ids, type=pa.string())]
                arrays += [pa.array(positive[:, j], mask=~finite[:, j], type=pa.bool_())
                           for j in range(len(markers))]
                writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
    return dict(raw_cells=size, measured_cells=measured_count,
                wholly_unmeasured_cells=missing_count,
                positive_counts=dict(zip(markers, counts.tolist())),
                output=str(destination))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", type=Path, default=FROZEN)
    parser.add_argument("--raw", type=Path, default=RAW)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("--output must be a new or empty directory")
    markers, gates = frozen_gates(args.frozen)
    args.output.mkdir(parents=True, exist_ok=True)
    donor_results = {}
    for donor in DONORS:
        source = args.raw / f"{donor}.h5ad"
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = args.output / f"{donor}-marker-calls.parquet"
        donor_results[donor] = stage_donor(source, destination, markers, gates,
                                           args.chunk_size)
        print(donor, donor_results[donor]["raw_cells"],
              donor_results[donor]["wholly_unmeasured_cells"], flush=True)
    report = dict(method="frozen_15_donor_log1p_raw_MFI_gates", frozen=str(args.frozen),
                  scope="marker positivity only; no phenotypes, scVI or app integration",
                  gate_values=dict(zip(markers, gates.tolist())), donors=donor_results)
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
