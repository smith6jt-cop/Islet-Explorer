#!/usr/bin/env python3
"""Stage and fit a joint 21-donor scVI model with the saved notebook recipe.

The 15 original cells come from the frozen phenotyped H5AD and are never
modified. Six additions enter the model only when they pass the *original*
notebook's QC rule: >=2 detected markers and total raw MFI at least the first
percentile of the original unfiltered cohort. Missing raw measurements are
excluded from training, but remain in the raw donor checkpoint.

The model uses raw MFI as its ``counts`` layer, Age/Gender categorical
covariates, 10 latent dimensions, two 128-wide layers, 0.1 dropout,
gene-batch dispersion and an NB likelihood. The notebook's prose mentions
``imageid`` batch correction, but its executable ``setup_anndata`` call omits
``batch_key``; this script follows the executable call. Its CLR/scale steps
affect visualization X, not the counts layer used by scVI.

    python scripts/retrain_joint_scvi.py prepare --output tmp/joint-scvi
    python scripts/retrain_joint_scvi.py train --output tmp/joint-scvi

This writes an offline derived candidate; it does not type cells, update the
app, or assert scientific acceptance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
UNFILTERED = ROOT.parent / "results/panc_CODEX.h5ad"
RAW = ROOT / "data/development/donor-additions-20260923-raw-v1/additions/raw"
DONORS = ("78", "6090", "6232", "6520", "6539", "6558")
CHUNK = 50_000


def text_values(values):
    return [value.decode() if isinstance(value, bytes) else str(value)
            for value in values]


def category_values(node, start, stop):
    if isinstance(node, h5py.Dataset):
        return text_values(node[start:stop])
    categories = text_values(node["categories"][:])
    return [categories[code] if code >= 0 else None
            for code in node["codes"][start:stop]]


def original_qc(path: Path, expected_kept: int):
    with h5py.File(path) as h5:
        matrix = h5["X"]
        totals = np.empty(matrix.shape[0], dtype=np.float64)
        detected = np.empty(matrix.shape[0], dtype=np.int8)
        for start in range(0, len(totals), CHUNK):
            end = min(start + CHUNK, len(totals))
            values = matrix[start:end]
            totals[start:end] = values.sum(axis=1)
            detected[start:end] = (values > 0).sum(axis=1)
    threshold = float(np.quantile(totals, 0.01))
    retained = (detected >= 2) & (totals >= threshold)
    if int(retained.sum()) != expected_kept:
        raise ValueError("Original notebook QC does not reproduce frozen cell count")
    return threshold, int(len(totals)), int(retained.sum())


def eligible(values: np.ndarray, threshold: float):
    finite = np.isfinite(values)
    complete = finite.all(axis=1)
    absent = ~finite.any(axis=1)
    if not np.all(complete | absent):
        raise ValueError("Partially measured addition needs separate review")
    if np.any(values[finite] < 0):
        raise ValueError("Negative raw MFI in addition")
    return complete & (values.sum(axis=1) >= threshold) & \
        ((values > 0).sum(axis=1) >= 2)


def count_addition(path: Path, markers: list[str], threshold: float):
    with h5py.File(path) as h5:
        source_markers = text_values(h5["var/marker"][:])
        if len(source_markers) != len(markers) or set(source_markers) != set(markers):
            raise ValueError(f"{path.name}: panel differs from original")
        count = 0
        for start in range(0, h5["X"].shape[0], CHUNK):
            values = h5["X"][start:start + CHUNK]
            count += int(eligible(values, threshold).sum())
        return count


def prepare(args):
    import anndata as ad
    import pandas as pd

    args.output.mkdir(parents=True, exist_ok=True)
    input_path = args.output / "training_input.h5ad"
    if input_path.exists():
        raise FileExistsError(input_path)
    with h5py.File(args.frozen) as h5:
        original_n = h5["X"].shape[0]
        markers = text_values(h5["raw/var/_index"][:])
        if markers != text_values(h5["var/_index"][:]):
            raise ValueError("Original raw and canonical marker orders differ")
    threshold, unfiltered_n, retained_n = original_qc(args.unfiltered, original_n)
    addition_counts = {donor: count_addition(args.raw / f"{donor}.h5ad", markers,
                                              threshold) for donor in DONORS}
    total = original_n + sum(addition_counts.values())
    values = np.empty((total, len(markers)), dtype=np.float32)
    ages = np.empty(total, dtype=np.float64)
    genders = []
    ids = []
    imageids = np.empty(total, dtype=np.int32)
    original = np.empty(total, dtype=bool)
    cursor = 0
    with h5py.File(args.frozen) as h5:
        for start in range(0, original_n, CHUNK):
            end = min(start + CHUNK, original_n)
            size = end - start
            values[cursor:cursor + size] = h5["raw/X"][start:end]
            ages[cursor:cursor + size] = h5["obs/Age"][start:end]
            im = h5["obs/imageid"][start:end].astype(np.int32)
            imageids[cursor:cursor + size] = im
            objects = text_values(h5["obs/Object ID"][start:end])
            ids.extend(f"{donor}:{obj}" for donor, obj in zip(im, objects))
            genders.extend(category_values(h5["obs/Gender"], start, end))
            original[cursor:cursor + size] = True
            cursor += size
    for donor in DONORS:
        path = args.raw / f"{donor}.h5ad"
        with h5py.File(path) as h5:
            source_markers = text_values(h5["var/marker"][:])
            order = [source_markers.index(marker) for marker in markers]
            for start in range(0, h5["X"].shape[0], CHUNK):
                end = min(start + CHUNK, h5["X"].shape[0])
                batch = h5["X"][start:end]
                use = eligible(batch, threshold)
                size = int(use.sum())
                if not size:
                    continue
                values[cursor:cursor + size] = batch[use][:, order]
                ages[cursor:cursor + size] = h5["obs/Age"][start:end][use]
                imageids[cursor:cursor + size] = int(donor)
                cells = text_values(h5["obs/cell_id"][start:end])
                ids.extend(cell for cell, keep in zip(cells, use) if keep)
                genders.extend(gender for gender, keep in
                               zip(category_values(h5["obs/Gender"], start, end), use)
                               if keep)
                original[cursor:cursor + size] = False
                cursor += size
        print(f"prepared donor {donor}: {addition_counts[donor]:,} cells", flush=True)
    if cursor != total or len(ids) != total or len(set(ids)) != total:
        raise ValueError("Training input cell census or IDs are inconsistent")
    if not np.all(np.isfinite(values)):
        raise ValueError("Training input contains missing measurements")
    obs = pd.DataFrame({"Age": ages, "Gender": genders, "imageid": imageids,
                        "original_donor": original}, index=pd.Index(ids, name="cell_id"))
    adata = ad.AnnData(X=values, obs=obs, var=pd.DataFrame(index=markers))
    adata.layers["counts"] = values
    adata.write_h5ad(input_path, compression="lzf")
    report = dict(method="SCVI_CODEX_v2.ipynb executable training settings",
                  source_original=str(args.frozen), source_unfiltered=str(args.unfiltered),
                  source_additions=str(args.raw), original_qc_percentile=0.01,
                  original_total_mfi_threshold=threshold,
                  original_unfiltered_cells=unfiltered_n,
                  original_retained_cells=retained_n,
                  addition_retained_cells=addition_counts, total_training_cells=total,
                  excluded_addition_cells_not_deleted_from_raw=True,
                  output=str(input_path), status="prepared_not_trained")
    (args.output / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"training input: {total:,} cells x {len(markers)} markers", flush=True)


def train(args):
    import anndata as ad
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    import scvi
    import torch

    input_path = args.output / "training_input.h5ad"
    model_path = args.output / "model"
    latent_path = args.output / "joint_scvi_latent.parquet"
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if model_path.exists() or latent_path.exists():
        raise FileExistsError("Model or latent output already exists")
    scvi.settings.seed = 0
    torch.set_num_threads(args.workers)
    adata = ad.read_h5ad(input_path)
    scvi.model.SCVI.setup_anndata(adata, layer="counts",
                                  categorical_covariate_keys=["Age", "Gender"],
                                  continuous_covariate_keys=None)
    model = scvi.model.SCVI(adata, n_latent=10, n_layers=2, n_hidden=128,
                            dropout_rate=0.1, dispersion="gene-batch",
                            gene_likelihood="nb")
    model.train(max_epochs=args.max_epochs, accelerator="gpu", devices=1,
                batch_size=args.batch_size, train_size=0.9,
                early_stopping=True, early_stopping_patience=45,
                check_val_every_n_epoch=1, plan_kwargs={"lr": 0.001},
                datasplitter_kwargs={"num_workers": args.workers,
                                     "pin_memory": True})
    model.save(model_path, overwrite=False, save_anndata=False)
    latent = model.get_latent_representation(adata, batch_size=args.batch_size)
    if latent.shape != (adata.n_obs, 10) or not np.all(np.isfinite(latent)):
        raise ValueError("Invalid joint scVI latent representation")
    fields = [pa.array(adata.obs_names.tolist(), type=pa.string())]
    fields += [pa.array(latent[:, j], type=pa.float32()) for j in range(10)]
    names = ["cell_id"] + [f"scVI_{j}" for j in range(10)]
    pq.write_table(pa.Table.from_arrays(fields, names=names), latent_path,
                   compression="zstd")
    history = getattr(model, "history", None)
    if history:
        pd.concat({name: frame for name, frame in history.items()}, axis=1).to_csv(
            args.output / "training_history.csv")
    report = dict(status="trained_derived_candidate_not_app_release",
                  model=str(model_path), latent=str(latent_path),
                  input_cells=adata.n_obs, latent_dimensions=latent.shape[1],
                  scvi_version=scvi.__version__, torch_version=torch.__version__,
                  max_epochs=args.max_epochs, batch_size=args.batch_size,
                  workers=args.workers, notebook_batch_key=None,
                  categorical_covariates=["Age", "Gender"])
    (args.output / "training.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "train"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--output", type=Path, required=True)
        if name == "prepare":
            cmd.add_argument("--frozen", type=Path, default=FROZEN)
            cmd.add_argument("--unfiltered", type=Path, default=UNFILTERED)
            cmd.add_argument("--raw", type=Path, default=RAW)
        else:
            cmd.add_argument("--max-epochs", type=int, default=400)
            cmd.add_argument("--batch-size", type=int, default=2048)
            cmd.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    else:
        train(args)


if __name__ == "__main__":
    main()
