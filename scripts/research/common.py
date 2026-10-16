"""Shared, explicit cohort and artifact contracts. No source discovery fallback."""
from __future__ import annotations
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from excluded_donors import assert_no_excluded

MARKERS = ['CD31','CD8a','CD3e','SMA','CD4','CD34','HLADR','PDPN','panCK','ECAD',
           'CD163','SST','ColIV','VIM','CD20','LGALS3','B3TUBB','GCG','KRT14','GAP43',
           'CD35','CHGA','PGP9.5','INS','CD44','CD45','NaKATPase','BCatenin','CD68','BActin','CK19']
RETIRED = ['Endocrine','Acinar','Endothelial','T cell','Immune','Unknown','ECAD+','Structural','APCs']
RETAINED = ['Alpha cell','Beta cell','Delta cell','B cell','CD4 Tcell','CD8a Tcell',
            'Macrophage','Neural','SMA+','Blood Vessel','Lymphatic','Ductal']
ENDOCRINE = ['Alpha cell','Beta cell','Delta cell']
SEED = 20260916

def require_fresh_output(path, sources=()):
    """Validate before creating anything; reject aliases and recursive staging."""
    path = Path(path)
    if path.is_symlink() or path.exists():
        raise FileExistsError(f'Fresh output root required: {path}')
    output = path.resolve()
    for source in sources:
        source = Path(source).resolve(strict=True)
        if output == source or output.is_relative_to(source) or source.is_relative_to(output):
            raise ValueError(f'Output must be separate from source: {source}')
    return output

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Replacing a file is safe even when a staged input was hardlinked.
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.'+path.name+'.', delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(json.dumps(value, indent=2, allow_nan=False, default=str)+'\n')
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)

def read_cohort(path, input_root):
    d = pd.read_csv(path, dtype=str).fillna('')
    # The historical 21-donor manifests name donor 78, removed 2026-10-15.
    assert_no_excluded(d.get('case_id', []), f'Cohort {Path(path).name}')
    required = {'case_id','donor_source','donor_status','cells','cell_slot','geometry','image','regional_root','metadata'}
    if not required.issubset(d) or d.case_id.duplicated().any():
        raise ValueError('Incomplete or duplicate cohort manifest')
    if len(d) != 21 or d.donor_status.value_counts().to_dict() != {'ND':7,'Aab+':7,'T1D':7}:
        raise ValueError('Expected exactly seven ND, seven Aab+, seven T1D donors')
    if (d.donor_source == 'HPAP').sum() != 1 or d.loc[d.donor_source == 'HPAP','case_id'].tolist() != ['78']:
        raise ValueError('HPAP identity must remain explicit')
    root = Path(input_root).resolve(strict=True)
    for col in ['cells','geometry','image','regional_root','metadata']:
        d[col] = d[col].map(lambda p: str((root / p).resolve(strict=True)))
    return d

def require_cohort(frame, donors, column='case_id', complete=True):
    if column not in frame:
        raise ValueError(f'Missing cohort key {column}')
    actual = set(frame[column].astype(str)); expected = set(map(str,donors))
    if actual - expected or (complete and expected - actual):
        raise ValueError(f'Mixed or incomplete cohort in {column}: {sorted(actual ^ expected)}')

def ratio(numerator, denominator):
    n, d = np.asarray(numerator, float), np.asarray(denominator, float)
    raw = np.divide(n, d, out=np.full(np.broadcast_shapes(n.shape,d.shape), np.nan),
                    where=np.isfinite(d) & (d > 0))
    return raw, np.log1p(raw)

def offline_report_root(input_root, report_root):
    """Analysis reports are explicit outputs outside the app's data root."""
    root, output = Path(input_root).resolve(), Path(report_root).resolve()
    if output == root or output.is_relative_to(root):
        raise ValueError('Research reports must be outside the application data root')
    output.mkdir(parents=True, exist_ok=True)
    return output
