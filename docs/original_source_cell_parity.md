# Original single-cell parity gate

`scripts/research/audit_original_source_cells.py` compares a future 21-donor
single-cell H5AD with the frozen
`single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` without
changing either file. It has not certified a 21-donor candidate. The
2026-10-03 protein expansion does not build a 21-donor single-cell H5AD: it
leaves the frozen file untouched (checked by SHA-256) and appends islet-level
rows, which `audit_original_parity.py` checks. This tool is kept for any future
single-cell rebuild. Run it only after such a candidate exists:

```bash
/home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  scripts/research/audit_original_source_cells.py \
  --candidate /path/to/21-donor-cells.h5ad \
  --report /path/outside/app/data/original_source_cell_parity.json
```

The checker aligns original cells by the pair `(obs/imageid, obs/Object ID)`;
the frozen `obs/_index` is a sequential row number, not a durable cell ID. It
requires every original cell exactly once and rejects added cells assigned to
an original donor. Candidate rows, markers and categorical label codes may be
reordered, and new donors or markers may be added.

With zero numerical tolerance, it checks frozen marker values in `X`, every
original layer including `layers/log`, `raw/X`, missingness, all original `obs`
source columns, phenotype labels, and the original donor-marker numeric
`uns/gates` table. Only the synthetic index and scVI registry columns are
excluded from protected `obs` comparison. Because a single-cell refit would
recompute them, `obsm/X_scVI`, `obsm/X_umap`, `obs/_scvi_batch` and
`obs/_scvi_labels` may change; every observed change or missing derived field
appears in `derived_differences` and its count in the CLI summary. A passing
source check does not accept the refit's biology, typing quality for the six
additions, or the islet-level app artifacts; those require separate review.

The implementation reads HDF5 fields directly in row chunks. It holds one
aligned original-cohort field at a time (about 627 MiB for a 2.65-million-row,
31-marker float64 matrix), rather than materializing the full multi-gigabyte
AnnData object. The report is created with exclusive file creation outside app
data so a prior report is not overwritten silently.
