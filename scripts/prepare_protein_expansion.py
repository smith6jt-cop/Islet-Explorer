#!/usr/bin/env python3
"""Prepare the added donors' measured protein inputs and a concrete rule review.

This does not assign addition phenotypes or publish app data. Owner-reviewed
rules and boundaries are required by the assembly step. Original data is read
only and its frozen phenotype is used strictly as a diagnostic comparator.
"""
from __future__ import annotations
import argparse
import html
import json
from pathlib import Path
import shutil

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from audit_phenotype_rule_replay import select_stratified_indices, agreement
from build_protein_embedding import sha256
from protein_phenotyping import RULES, read_rules, classify, qc_mask, rescale_proteins
from stage_original_gate_calls import frozen_gates, DONORS


def read_column(h5, name, rows=slice(None)):
    node = h5['obs'][name]
    if isinstance(node, h5py.Group):
        cats = node['categories'].asstr()[:] if node['categories'].dtype.kind in 'OSU' else node['categories'][:]
        codes = node['codes'][rows]
        return np.array([cats[c] if c >= 0 else None for c in codes])
    return node.asstr()[rows] if node.dtype.kind in 'OSU' else node[rows]


def sample_matrix(node, rows):
    out = np.empty((len(rows), node.shape[1]), dtype=node.dtype)
    for start in range(0, node.shape[0], 50000):
        use = np.flatnonzero((rows >= start) & (rows < start + 50000))
        if len(use):out[use] = node[start:start+50000][rows[use] - start]
    return out


def review_originals(frozen, out, rules):
    with h5py.File(frozen) as h5:
        phenotypes = h5['obs/phenotype/categories'].asstr()[:].tolist()
        donors = h5['obs/donor_id/categories'].asstr()[:].tolist()
        pc, dc = h5['obs/phenotype/codes'][:], h5['obs/donor_id/codes'][:]
        rows, weights, strata = select_stratified_indices(pc, dc, phenotypes, donors, 30, 42)
        markers = h5['var/_index'].asstr()[:].tolist()
        values = sample_matrix(h5['X'], rows)
        labels, reasons = classify(values, markers, rules)
        truth = np.asarray(phenotypes)[pc[rows]]
        donor = np.asarray(donors)[dc[rows]]
        frame = pd.DataFrame(dict(donor=donor, saved=truth, proposed=labels, reason=reasons,
            weight=weights, Object_ID=read_column(h5, 'Object ID', rows),
            x_um=read_column(h5, 'X_centroid', rows), y_um=read_column(h5, 'Y_centroid', rows)))
    frame.to_csv(out / 'original_sample.csv', index=False)
    pd.crosstab(frame.saved, frame.proposed).to_csv(out / 'confusion.csv')
    records=[]
    for column in ('donor','saved'):
        for label, group in frame.groupby(column):
            records.append(dict(grouping=column, group=label, **agreement(group.saved.to_numpy(),
                group.proposed.to_numpy(), group.weight.to_numpy())))
    pd.DataFrame(records).to_csv(out / 'agreement_by_donor_and_label.csv', index=False)
    marker_frame = pd.DataFrame(values, columns=markers)
    marker_frame['saved'] = truth
    marker_frame.groupby('saved').median().to_csv(out / 'saved_label_marker_medians.csv')
    marker_frame['proposed'] = labels
    marker_frame.groupby('proposed')[markers].median().to_csv(out / 'proposed_label_marker_medians.csv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    table = marker_frame.groupby('saved')[markers].median()
    fig, ax = plt.subplots(figsize=(15, 8))
    im=ax.imshow(table.to_numpy(), vmin=0, vmax=1, cmap='viridis', aspect='auto')
    ax.set_xticks(range(len(markers)), markers, rotation=90)
    ax.set_yticks(range(len(table)), table.index)
    ax.set_title('Frozen labels: median gate-rescaled protein signal (stratified sample)')
    fig.colorbar(im, ax=ax, label='Gate-rescaled signal; positive threshold 0.5')
    fig.tight_layout();fig.savefig(out / 'original_marker_evidence.png', dpi=150);plt.close(fig)
    result=dict(sample_cells=len(rows), sampling='30 per saved label × donor, seed 42',
        exact_reproduction_required=False, **agreement(truth, labels, weights),
        proposed_label_counts=frame.proposed.value_counts().to_dict(),
        reason_counts=frame.reason.value_counts().to_dict(),
        unsupported_historical_labels=sorted(set(truth) - set(labels)))
    (out / 'rule_diagnostics.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def prepare(frozen, raw_root, output, report_root):
    output, report_root = Path(output).resolve(), Path(report_root).resolve()
    if output.exists() or report_root.exists():
        raise FileExistsError('Use fresh output and report directories')
    if report_root.is_relative_to(output):
        raise ValueError('Review evidence must be outside the data bundle')
    for source in (Path(frozen).resolve(),Path(raw_root).resolve()):
        if output == source or output.is_relative_to(source) or source.is_relative_to(output):
            raise ValueError('Output must be separate from inputs')
    markers, gates = frozen_gates(frozen)
    rules = read_rules()
    output.mkdir(parents=True);report_root.mkdir(parents=True)
    raw_digests={cid:sha256(Path(raw_root)/f'{cid}.h5ad') for cid in DONORS}
    frozen_digest=sha256(frozen)
    records=[]
    for cid in DONORS:
        path=Path(raw_root)/f'{cid}.h5ad'
        with h5py.File(path) as h5:
            source_markers=h5['var/marker'].asstr()[:].tolist()
            if set(source_markers)!=set(markers):raise ValueError(f'{cid}: marker panel mismatch')
            raw=h5['X'][:][:,[source_markers.index(m) for m in markers]]
            keep=qc_mask(raw)
            ids=read_column(h5,'cell_id')
            if len(set(ids)) != len(ids):raise ValueError(f'{cid}: duplicate source identities')
            columns=['Object ID','Parent','X_centroid','Y_centroid','Cell Area','Nucleus Area','Age','Gender','imageid','donor_id','donor_source']
            obs=pd.DataFrame({c:read_column(h5,c) for c in columns if c in h5['obs']},index=ids)
        measured=np.isfinite(raw).all(axis=1)
        exclusion=np.where(~measured,'wholly_unmeasured',np.where(~keep,'below_original_qc','included'))
        pd.DataFrame(dict(cell_id=ids,case_id=cid,analysis_included=keep,reason=exclusion)).to_parquet(
            output / f'{cid}-census.parquet',index=False)
        scaled,logged,bounds=rescale_proteins(raw[keep], gates)
        cells=ad.AnnData(scaled,obs=obs.iloc[np.flatnonzero(keep)].copy(),var=pd.DataFrame(index=markers))
        cells.raw=ad.AnnData(raw[keep],obs=cells.obs.copy(),var=cells.var.copy())
        cells.layers['log']=logged
        cells.uns['protein_preparation']=dict(source_sha256=raw_digests[cid],
            source=str(path.resolve()),scale='per-donor gate-rescaled log1p MFI',**bounds)
        cells.uns['gates']=pd.DataFrame({'markers':markers,cid:gates})
        cells.write_h5ad(output/f'{cid}.h5ad',compression='gzip')
        records.append(dict(case_id=cid,raw_cells=len(raw),measured_cells=int(measured.sum()),
            qc_retained=int(keep.sum()),wholly_unmeasured=int((~measured).sum()),
            measured_qc_excluded=int((measured & ~keep).sum()),
            prepared_sha256=sha256(output/f'{cid}.h5ad'),
            markers_gate_outside_observed_range=[m for m,flag in zip(markers,bounds['gate_outside_observed_range']) if flag]))
        print(f'{cid}: {keep.sum():,} QC-retained; no phenotypes assigned',flush=True)
    diagnostics=review_originals(frozen,report_root,rules)
    pd.DataFrame(records).to_csv(report_root/'donor_census.csv',index=False)
    shutil.copy2(RULES,report_root/RULES.name)
    manifest=dict(status='scaled_untyped_waiting_for_rule_and_boundary_review',
        original_source=str(Path(frozen).resolve()),original_source_sha256=frozen_digest,
        raw_source_sha256=raw_digests,rules_sha256=sha256(RULES),donors=records,
        app_ready=False,phenotypes_assigned=False,
        next_step=f'Review rules and {len(DONORS)} donor boundary decisions before assembly')
    (output/'preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
    rules_frame=pd.read_csv(RULES).fillna('')
    (report_root/'index.html').write_text('''<!doctype html><html lang="en"><meta charset="utf-8">
<title>Protein expansion review</title><style>body{font:16px system-ui;max-width:1400px;margin:2rem auto;padding:1rem}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:6px}img{max-width:100%}</style>
<h1>Protein expansion: rules and boundary review</h1>
<p>This proposal is not the historical phenotype method. Original labels and measurements remain unchanged. No addition phenotypes have been assigned.</p>
<p>Positive means gate-rescaled signal ≥0.5. Competing top-level lineages become Unknown; competing child types retain their parent. Co-positive hormones retain Endocrine. Acinar, ECAD+ and Structural remain unsupported instead of receiving invented negative-only identities.</p>
<h2>Proposed rules</h2>'''+rules_frame.to_html(index=False,escape=True)+
        '<h2>Original-donor diagnostic</h2><p>Sample agreement: '+f"{diagnostics['sample_agreement']:.1%}"+
        '. This is a stratified diagnostic, not an acceptance threshold or a population accuracy estimate.</p>'+
        '<p><a href="confusion.csv">Confusion table</a> · <a href="agreement_by_donor_and_label.csv">Per-donor and label results</a> · <a href="original_sample.csv">Representative cell IDs and coordinates</a></p>'+
        '<img src="original_marker_evidence.png" alt="Median marker evidence by frozen phenotype">'+
        '<h2>Added-donor QC census</h2>'+pd.DataFrame(records).drop(columns=['prepared_sha256']).to_html(index=False,escape=True)+
        '<h2>Boundary review</h2><p>Each donor needs an explicit owner boundary decision in the acceptance record. Use the current native QuPath project and the existing gallery linked in the execution record.</p></html>')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frozen',type=Path,required=True);p.add_argument('--raw-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--report-root',type=Path,required=True)
    a=p.parse_args();prepare(a.frozen,a.raw_root,a.output,a.report_root)
