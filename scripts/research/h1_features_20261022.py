#!/usr/bin/env python3
"""Read-only H1 feature extraction from the accepted uniform v10 app root.

No models are fitted here. All identities use persistent object_id/case_id joins,
never row-order joins. Geometry is converted from image pixels to micrometers.
Primary peri counts retain the app's overlapping object-associated bands; an
exclusive nearest-Islet/Cluster allocation is a prespecified sensitivity.

Example:
  /home/smith6jt/miniconda3/envs/scvi-env/bin/python \
    scripts/research/h1_features_20261022.py \
    --data-root data/development/uniform-20donor-20261022-v10 \
    --output-root output/h1_20261022 --workers 2
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import hashlib
import io
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from shapely.affinity import scale
from shapely.geometry import shape

PIXEL_UM = 0.3774
LABELS = ["Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative", "Endocrine",
          "CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune", "Exocrine"]
IMMUNE = {"CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune"}
ENDOCRINE = set(LABELS[:5])
OBJECT_CLASSES = {"Islet", "EndocrineCluster"}


def token(label):
    return label.replace(" ", "_").replace("+", "plus")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite_json(x):
    if isinstance(x, dict):
        return {str(k): finite_json(v) for k, v in x.items()}
    if isinstance(x, (tuple, list, np.ndarray)):
        return [finite_json(v) for v in x]
    if isinstance(x, (np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


def write_json(path, x):
    Path(path).write_text(json.dumps(finite_json(x), indent=2, allow_nan=False) + "\n")


def density(n, area_um2):
    return n * 1e6 / area_um2 if area_um2 > 0 else np.nan


def counts(labels, ids):
    return np.bincount(labels[ids], minlength=len(LABELS))


def union_local(tree, geoms, target):
    ix = tree.query(target, predicate="intersects")
    return shapely.union_all(geoms[ix]) if len(ix) else shapely.GeometryCollection()


def geometry(root, case):
    path = root / "gson" / f"{case:04d}.geojson.gz"
    with gzip.open(path, "rt") as f:
        obj = json.load(f)
    cores, expanded, tissue, invalid = {}, {}, [], 0
    for feature in obj["features"]:
        p = feature.get("properties") or {}
        cl = (p.get("classification") or {}).get("name")
        if cl not in {"Tissue", "Islet", "IsletExpanded"}:
            continue
        if cl != "Tissue" and p.get("object_class") not in OBJECT_CLASSES:
            raise ValueError(f"{case}: unexpected analysis class {p}")
        g = shape(feature["geometry"])
        if not g.is_valid:
            invalid += 1
            g = shapely.make_valid(g)
        g = scale(g, xfact=PIXEL_UM, yfact=PIXEL_UM, origin=(0, 0))
        # Floating-point rescaling can turn a touching ring into an invalid
        # ring even when the saved pixel polygon was valid (observed on6090).
        # Repair only the in-memory derived geometry, then verify cell counts.
        if not g.is_valid:
            invalid += 1
            g = shapely.make_valid(g)
        if cl == "Tissue":
            tissue.append(g)
        elif cl == "Islet":
            cores[p["name"]] = (p["object_class"], g)
        else:
            expanded[p["name"].removesuffix("_exp20um")] = g
    del obj
    if not tissue or set(cores) != set(expanded):
        raise ValueError(f"{case}: missing tissue or core/expanded pairing")
    return cores, expanded, shapely.union_all(tissue), invalid, path


def exclusive_peri_owners(points, core_geoms, names, peri_ids, eligible_objects=None):
    """Nearest boundary among actual band memberships; stable name breaks ties.

    The many-to-many bands remain the primary H1 exposure. Only additive
    censuses and the existing exclusive-peri sensitivity use this ownership.
    """
    owner = np.full(len(points), -1, dtype=np.int32)
    best = np.full(len(points), np.inf)
    eligible_objects = np.ones(len(names), bool) if eligible_objects is None else eligible_objects
    for i in sorted(range(len(names)), key=lambda j: names[j]):
        if not eligible_objects[i]:
            continue
        ids = peri_ids[i]
        dist = shapely.distance(points[ids], core_geoms[i].boundary)
        wins = dist < best[ids] - 1e-9
        chosen = ids[wins]
        owner[chosen], best[chosen] = i, dist[wins]
    return owner


def extract_donor(root_string, output_string, case, obs, metadata, verify_drilldowns, eligibility_root=None):
    start = time.monotonic()
    root, out = Path(root_string), Path(output_string)
    cores, expanded, tissue, invalid, geom_path = geometry(root, case)
    names = list(cores)
    core_geoms = np.asarray([cores[n][1] for n in names], dtype=object)
    expanded_geoms = np.asarray([expanded[n] for n in names], dtype=object)
    bands = np.asarray([expanded[n].difference(cores[n][1]) for n in names], dtype=object)
    ctree, etree = STRtree(core_geoms), STRtree(expanded_geoms)
    name_ix = {n: i for i, n in enumerate(names)}
    paired = set(obs.base_islet_id.astype(str))
    cell_path = root / "cell_types" / f"{case}.parquet"
    cells = pd.read_parquet(cell_path)
    if cells.object_id.duplicated().any() or set(cells.case_id.astype(str)) != {str(case)}:
        raise ValueError(f"{case}: duplicate cell identity or incorrect donor")
    if not set(cells.phenotype).issubset(LABELS):
        raise ValueError(f"{case}: unknown phenotype")
    xy = cells[["x_um", "y_um"]].to_numpy(float)
    if not np.isfinite(xy).all():
        raise ValueError(f"{case}: invalid cell coordinates")
    points = shapely.points(xy)
    ptree = STRtree(points)
    lab = pd.Categorical(cells.phenotype, categories=LABELS).codes.astype(np.int16)
    immune = cells.phenotype.isin(IMMUNE).to_numpy()
    in_core = cells.cell_region.eq("core").to_numpy()
    # Persistent identity is authoritative for source core assignment. Geometry
    # rederivation is checked below, and never silently replaces these counts.
    memberships = cells.islet_name.fillna("").astype(str).to_numpy()
    core_owner = np.full(len(cells), -1, dtype=np.int32)
    for n, i in name_ix.items():
        core_owner[in_core & (memberships == n)] = i
    if np.any(in_core & (core_owner < 0)):
        raise ValueError(f"{case}: core cells have no geometric owner")

    cell_keep = np.ones(len(cells), bool)
    object_keep = np.ones(len(names), bool)
    edge_objects, edge_paths = None, []
    if eligibility_root:
        edge_paths = [Path(eligibility_root) / 'cells' / f'{case}.parquet',
                      Path(eligibility_root) / 'objects' / f'{case}.parquet',
                      Path(eligibility_root) / 'edge_manifest.json']
        flags = pd.read_parquet(edge_paths[0]).set_index('object_id')
        if not flags.index.is_unique or not cells.object_id.isin(flags.index).all():
            raise ValueError(f'{case}: incomplete/duplicate edge-cell join')
        cell_keep = flags.loc[cells.object_id, 'edge_analysis_eligible'].to_numpy(bool)
        edge_objects = pd.read_parquet(edge_paths[1]).set_index('object_name')
        if not edge_objects.index.is_unique or not set(names).issubset(edge_objects.index):
            raise ValueError(f'{case}: incomplete/duplicate edge-object join')
        object_keep = edge_objects.loc[names, 'core_complete'].to_numpy(bool)

    peri_ids = []
    reuse = np.zeros(len(cells), dtype=np.int16)
    core_geom_mismatch = []
    for i, g in enumerate(core_geoms):
        ci_geom = ptree.query(g, predicate="covers")
        ci = np.flatnonzero(core_owner == i)
        if not np.array_equal(np.sort(ci_geom), ci):
            core_geom_mismatch.append({"name": names[i], "geometry": len(ci_geom), "identity": len(ci),
                                       "symmetric_difference": len(np.setxor1d(ci_geom, ci))})
        ids = ptree.query(bands[i], predicate="covers")
        ids = ids[~in_core[ids]]
        peri_ids.append(ids)
    raw_peri_ids = peri_ids
    peri_ids = [ids[cell_keep[ids]] for ids in raw_peri_ids]
    for i, ids in enumerate(peri_ids):
        if object_keep[i]:
            reuse[ids] += 1
    union_peri = np.flatnonzero(reuse > 0)
    # Nearest boundary among ALL Islets/Clusters, including unpaired objects.
    # Singles never compete for an assignment. Only actual peri-band members
    # are assigned, so tissue clipping/other-core exclusions remain respected.
    nearest_owner = exclusive_peri_owners(points, core_geoms, names, peri_ids, object_keep)
    # Explicit check rather than assuming the closest object supplied a band.
    covered_exclusive = np.zeros(len(cells), bool)
    for i, ids in enumerate(peri_ids):
        covered_exclusive[ids[nearest_owner[ids] == i]] = True
    if not covered_exclusive[union_peri].all():
        raise ValueError(f"{case}: nearest object does not cover one of its allocated peri cells")

    all_exclusion = shapely.union_all(expanded_geoms)
    distal_region = tissue.difference(all_exclusion)
    distal_ids = ptree.query(distal_region, predicate="covers")
    distal_ids = distal_ids[cell_keep[distal_ids]]
    geometrically_in_tissue = shapely.covers(tissue, points)
    inside_tissue = geometrically_in_tissue & cell_keep
    tissue_area = float(tissue.area)
    tissue_count = int(immune[inside_tissue].sum())
    distal_immune = int(immune[distal_ids].sum())
    rows, discrepancies, drilldown_hashes = [], [], {}
    for rec in obs.to_dict("records"):
        name = str(rec["base_islet_id"])
        i = name_ix[name]
        core, band = core_geoms[i], bands[i]
        raw_ci, raw_pi = np.flatnonzero(core_owner == i), raw_peri_ids[i]
        ci, pi = raw_ci[cell_keep[raw_ci]], peri_ids[i]
        exc = pi[nearest_owner[pi] == i]
        cc, pc = counts(lab, ci), counts(lab, pi)
        raw_cc, raw_pc = counts(lab, raw_ci), counts(lab, raw_pi)
        for j, label in enumerate(LABELS):
            expected_core = int(rec["count_" + label])
            expected_peri = int(rec["peri_count_" + token(label)])
            if (int(raw_cc[j]), int(raw_pc[j])) != (expected_core, expected_peri):
                discrepancies.append({"name": name, "phenotype": label, "core_derived": int(raw_cc[j]),
                                      "core_app": expected_core, "peri_derived": int(raw_pc[j]), "peri_app": expected_peri})
        n_core_immune, n_peri_immune = int(immune[ci].sum()), int(immune[pi].sum())
        n_exc_immune = int(immune[exc].sum())
        n_near, n_exc_near = n_core_immune + n_peri_immune, n_core_immune + n_exc_immune
        row = dict(islet_id=rec["islet_id"], case_id=case, donor_status=str(rec["donor_status"]),
                   object_class=str(rec["object_class"]), base_islet_id=name,
                   core_area_um2=float(core.area), peri_area_um2=float(band.area),
                   core_perimeter_um=float(core.length), centroid_x_um=float(core.centroid.x),
                   centroid_y_um=float(core.centroid.y), n_cells_core=len(ci), n_cells_peri=len(pi),
                   core_tissue_fraction=float(core.intersection(tissue).area / core.area),
                   core_distance_to_tissue_edge_um=float(core.distance(tissue.boundary)),
                   n_beta_core=int(cc[0]), n_endocrine_core=int(sum(cc[j] for j, l in enumerate(LABELS) if l in ENDOCRINE)),
                   n_immune_core=n_core_immune, n_immune_peri=n_peri_immune, n_immune_nearby=n_near,
                   immune_localization_fraction=n_core_immune / n_near if n_near else np.nan,
                   immune_exposed=n_near > 0, excluded_cells_core=int(rec["excluded_cells_core"]),
                   n_cells_peri_exclusive=len(exc), n_immune_peri_exclusive=n_exc_immune,
                   n_immune_nearby_exclusive=n_exc_near,
                   immune_localization_fraction_exclusive=n_core_immune / n_exc_near if n_exc_near else np.nan,
                   immune_exposed_exclusive=n_exc_near > 0,
                   peri_reused_cell_count=int((reuse[pi] > 1).sum()),
                   peri_mean_membership_multiplicity=float(reuse[pi].mean()) if len(pi) else np.nan,
                   tissue_area_mm2=tissue_area / 1e6, tissue_cell_count=int(inside_tissue.sum()),
                   tissue_immune_count=tissue_count, tissue_immune_density_mm2=density(tissue_count, tissue_area),
                   analysed_cell_count_all=int(cell_keep.sum()), analysed_immune_count_all=int(immune[cell_keep].sum()),
                   analysed_cells_outside_tissue=int((cell_keep & ~geometrically_in_tissue).sum()),
                   analysed_immune_outside_tissue=int(immune[cell_keep & ~geometrically_in_tissue].sum()),
                   distal_tissue_area_um2=float(distal_region.area), distal_tissue_cell_count=len(distal_ids),
                   distal_tissue_immune_count=distal_immune,
                   distal_tissue_immune_density_mm2=density(distal_immune, distal_region.area))
        row.update({f"n_core_{token(l)}": int(cc[j]) for j, l in enumerate(LABELS)})
        row.update({f"n_peri_{token(l)}": int(pc[j]) for j, l in enumerate(LABELS)})
        if edge_objects is not None:
            er = edge_objects.loc[name]
            row.update(core_complete=bool(er.core_complete), edge_exclusion_reason=er.exclusion_reason,
                       edge_cells_removed_core=len(raw_ci)-len(ci), edge_cells_removed_peri=len(raw_pi)-len(pi),
                       source_n_cells_core=len(raw_ci), source_n_cells_peri=len(raw_pi))
            for radius in (20, 100, 250):
                for prefix in ('neighborhood_complete_', 'neighborhood_requested_area_', 'neighborhood_observed_area_'):
                    key = prefix + str(radius) + ('' if prefix == 'neighborhood_complete_' else '_um2')
                    row[key] = er[key]
        for lo, hi in [(20, 100), (100, 250)]:
            outer = core.buffer(hi)
            full_ring = outer.difference(core.buffer(lo))
            ring = full_ring.intersection(tissue)
            tissue_fraction = float(ring.area / full_ring.area) if full_ring.area else np.nan
            # Subtract every object's immediate core+20um peri footprint. The
            # local union avoids repeatedly overlaying the entire image.
            ring = ring.difference(union_local(etree, expanded_geoms, outer))
            ids = ptree.query(ring, predicate="covers") if not ring.is_empty else np.array([], dtype=int)
            ids = ids[cell_keep[ids]]
            pref = f"field{lo}_{hi}"
            row[pref + "_area_um2"] = float(ring.area)
            row[pref + "_tissue_fraction"] = tissue_fraction
            row[pref + "_background_fraction"] = float(ring.area / full_ring.area) if full_ring.area else np.nan
            row["n_cells_" + pref] = len(ids)
            row["n_immune_" + pref] = int(immune[ids].sum())
            row["immune_density_" + pref + "_mm2"] = density(immune[ids].sum(), ring.area)
            row[pref + "_assessable"] = bool(ring.area > 0)
            if edge_objects is not None and not row[f'neighborhood_complete_{hi}']:
                row['immune_density_' + pref + '_mm2'] = np.nan
                row[pref + '_assessable'] = False
            if np.any(in_core[ids]) or np.any(reuse[ids] > 0):
                raise ValueError(f"{case}/{name}: local field includes core/immediate-peri cells")
        row.update(metadata)
        if verify_drilldowns:
            path = root / "cells" / f"{rec['islet_id']}.csv"
            content = path.read_bytes()
            drilldown_hashes[str(path)] = hashlib.sha256(content).hexdigest()
            drill = pd.read_csv(io.BytesIO(content), usecols=["phenotype", "cell_region"])
            for region, expected in [("core", raw_cc), ("peri", raw_pc)]:
                actual = drill.loc[drill.cell_region.eq(region), "phenotype"].value_counts().reindex(LABELS, fill_value=0).to_numpy()
                if not np.array_equal(actual, expected):
                    raise ValueError(f"{case}/{name}: drilldown {region} count mismatch")
        rows.append(row)
    if discrepancies:
        write_json(out / f"FAILED_count_discrepancies_{case}.json", discrepancies)
        raise ValueError(f"{case}: geometry-derived peri/source core counts differ from served H5AD; see discrepancy report")
    features = pd.DataFrame(rows)
    features.to_parquet(out / "features_by_donor" / f"{case}.parquet", index=False)
    u = cells.iloc[union_peri][["object_id", "case_id", "phenotype", "x_um", "y_um"]].copy()
    u["membership_multiplicity"] = reuse[union_peri]
    u["exclusive_owner"] = np.asarray(names, dtype=object)[nearest_owner[union_peri]]
    u["exclusive_owner_in_paired_table"] = u.exclusive_owner.isin(paired)
    u.to_parquet(out / "peri_cells_by_donor" / f"{case}.parquet", index=False)
    if eligibility_root:
        primary = cells[['object_id', 'case_id', 'phenotype', 'x_um', 'y_um']].copy()
        primary['edge_analysis_eligible'] = cell_keep
        primary['primary_compartment'] = 'tissue'
        primary['primary_object'] = ''
        primary.loc[union_peri, 'primary_compartment'] = 'peri'
        primary.loc[union_peri, 'primary_object'] = np.asarray(names)[nearest_owner[union_peri]]
        for i, name in enumerate(names):
            ids = np.flatnonzero((core_owner == i) & cell_keep)
            primary.loc[ids, 'primary_compartment'] = 'core' if object_keep[i] else 'excluded_object'
            primary.loc[ids, 'primary_object'] = name
        primary.loc[~cell_keep, 'primary_compartment'] = 'excluded_cell'
        primary.to_parquet(out / 'analysis_cells_by_donor' / f'{case}.parquet', index=False)
        relations = []
        for i, name in enumerate(names):
            for region, ids in [('core', np.flatnonzero((core_owner == i) & cell_keep)), ('peri', peri_ids[i])]:
                r = cells.iloc[ids][['object_id', 'case_id', 'phenotype']].copy()
                r['islet_id'], r['cell_region'], r['core_complete'] = f'{case}_{name}', region, bool(object_keep[i])
                r['exclusive_owner'] = (nearest_owner[ids] == i) if region == 'peri' else object_keep[i]
                relations.append(r)
        relation = pd.concat(relations, ignore_index=True)
        if relation.duplicated(['case_id', 'object_id', 'islet_id']).any():
            raise ValueError('Duplicate physical cell within an object core/peri union')
        relation.to_parquet(out / 'object_memberships_by_donor' / f'{case}.parquet', index=False)
    audit = dict(case_id=case, donor_status=str(obs.donor_status.iloc[0]),
                 all_geometric_objects=len(names), paired_objects=len(rows), analysed_cells_all=len(cells),
                 tissue_cells=int(inside_tissue.sum()), analysed_immune_count_all=int(immune[cell_keep].sum()),
                 cells_outside_tissue_geometry=int((cell_keep & ~geometrically_in_tissue).sum()),
                 immune_outside_tissue_geometry=int(immune[cell_keep & ~geometrically_in_tissue].sum()),
                 edge_cells_excluded=int((~cell_keep).sum()), edge_objects_excluded=int((~object_keep).sum()),
                 eligible_unique_cells=int(cell_keep.sum()),
                 n_unique_peri_cells=len(union_peri), n_peri_memberships_all_objects=int(reuse.sum()),
                 n_unique_reused_peri_cells=int((reuse > 1).sum()), maximum_peri_reuse=int(reuse.max()),
                 core_geometry_identity_discrepancies=core_geom_mismatch,
                 count_discrepancies=discrepancies, repaired_geometry_features=invalid,
                 tissue_area_mm2=tissue_area / 1e6, tissue_immune_count=tissue_count,
                 tissue_immune_density_mm2=density(tissue_count, tissue_area),
                 distal_tissue_area_um2=float(distal_region.area), distal_tissue_immune_count=distal_immune,
                 distal_tissue_immune_density_mm2=density(distal_immune, distal_region.area),
                 elapsed_seconds=round(time.monotonic() - start, 2),
                 source_sha256={str(p): sha256(p) for p in [geom_path, cell_path, *edge_paths]},
                 drilldown_sha256=drilldown_hashes)
    write_json(out / "donor_audits" / f"{case}.json", audit)
    print(f"case {case}: {len(rows)} paired objects; all source/H5AD/drilldown counts agree; "
          f"{audit['n_unique_reused_peri_cells']} shared peri cells; {audit['elapsed_seconds']}s", flush=True)
    return audit


def verified_checkpoint(root, out, case, obs, metadata, require_drilldowns):
    """Reuse only source-hash-verified checkpoints whose counts match live obs."""
    ap, fp = out / "donor_audits" / f"{case}.json", out / "features_by_donor" / f"{case}.parquet"
    if not ap.exists() or not fp.exists():
        return None
    audit = json.loads(ap.read_text())
    sources = dict(audit["source_sha256"])
    if require_drilldowns:
        if len(audit.get("drilldown_sha256", {})) != len(obs):
            return None
        sources.update(audit["drilldown_sha256"])
    if any(sha256(p) != expected for p, expected in sources.items()):
        raise ValueError(f"{case}: checkpoint source hash changed")
    d = pd.read_parquet(fp).set_index("islet_id").sort_index()
    reference = obs.set_index("islet_id").sort_index()
    if not d.index.equals(reference.index):
        raise ValueError(f"{case}: checkpoint object keys changed")
    for label in LABELS:
        for dc, rc in [(f"n_core_{token(label)}", "count_" + label),
                       (f"n_peri_{token(label)}", "peri_count_" + token(label))]:
            if not np.array_equal(d[dc].to_numpy(), reference[rc].to_numpy()):
                raise ValueError(f"{case}: checkpoint phenotype counts changed")
    for column, expected in metadata.items():
        if pd.isna(expected):
            good = d[column].isna().all()
        else:
            good = d[column].eq(expected).all()
        if not good:
            raise ValueError(f"{case}: checkpoint donor metadata changed: {column}")
    print(f"case {case}: reused hash/count/metadata-verified checkpoint", flush=True)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--cases", type=int, nargs="+", help="Explicit pilot subset; final run must omit this")
    parser.add_argument("--skip-drilldown-verification", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Reuse only source-hash/count/metadata-verified donor checkpoints")
    parser.add_argument("--checkpoint-only", action="store_true", help="Write/validate donor checkpoints without aggregate tables")
    parser.add_argument('--eligibility-root', type=Path, help='Analysis-only cell/object edge audit; preserves baseline count verification')
    args = parser.parse_args()
    root, out = args.data_root.resolve(strict=True), args.output_root.resolve()
    if out == root or out.is_relative_to(root) or root.is_relative_to(out):
        raise ValueError("Outputs must be outside all app/source data")
    for sub in ["features_by_donor", "peri_cells_by_donor", "donor_audits"]:
        (out / sub).mkdir(parents=True, exist_ok=True)
    if args.eligibility_root:
        for sub in ('analysis_cells_by_donor', 'object_memberships_by_donor'):
            (out / sub).mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(root / "islet_explorer.h5ad", backed="r")
    obs = a.obs.copy()
    a.file.close()
    obs["case_id"] = obs.imageid.astype(int)
    donors = pd.read_csv(root / "donor_metadata.csv")
    census = pd.read_csv(root / "census.csv")
    if 78 in set(obs.case_id) or 78 in set(donors.case_id) or 78 in set(census.case_id):
        raise ValueError("Excluded donor 78 present")
    if len(donors) != 20 or set(obs.case_id) != set(donors.case_id) or obs.islet_id.duplicated().any():
        raise ValueError("Not the accepted complete 20-donor cohort")
    groups = obs.drop_duplicates("case_id").donor_status.value_counts().to_dict()
    if groups != {"ND": 7, "Aab+": 7, "T1D": 6}:
        raise ValueError(f"Unexpected donor groups: {groups}")
    cases = sorted(set(args.cases) if args.cases else set(donors.case_id))
    if not set(cases).issubset(set(donors.case_id)):
        raise ValueError("Requested case is absent/excluded")
    source_paths = [root / f for f in ["islet_explorer.h5ad", "donor_metadata.csv", "census.csv",
                                      "uniform_release.json", "islet_spatial_lookup.csv", "image_sources.json"]]
    source_hashes = {str(p): sha256(p) for p in source_paths}
    audits = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs = []
        for case in cases:
            md = donors.loc[donors.case_id.eq(case)].iloc[0].to_dict()
            md.pop("case_id")
            if args.resume and not args.eligibility_root:
                audit = verified_checkpoint(root, out, case, obs.loc[obs.case_id.eq(case)], md,
                                            not args.skip_drilldown_verification)
                if audit is not None:
                    audits.append(audit)
                    continue
            fs.append(pool.submit(extract_donor, str(root), str(out), case, obs.loc[obs.case_id.eq(case)],
                                  md, not args.skip_drilldown_verification, args.eligibility_root))
        for future in as_completed(fs):
            audits.append(future.result())
    if args.checkpoint_only:
        print(f"Checkpoint extraction complete for {len(cases)} donors; aggregate tables not written.", flush=True)
        return
    tables = [pd.read_parquet(out / "features_by_donor" / f"{case}.parquet") for case in cases]
    result = pd.concat(tables, ignore_index=True).sort_values(["case_id", "islet_id"]).reset_index(drop=True)
    if len(result) != int(obs.case_id.isin(cases).sum()) or result.islet_id.duplicated().any():
        raise ValueError("Incomplete/duplicated object feature table")
    result.to_csv(out / "h1_object_features.csv", index=False)
    result.to_parquet(out / "h1_object_features.parquet", index=False)
    donor_table = pd.DataFrame([{k: v for k, v in x.items() if k not in
                                {"source_sha256", "drilldown_sha256", "core_geometry_identity_discrepancies", "count_discrepancies"}}
                               for x in audits]).sort_values("case_id")
    donor_table = donor_table.merge(donors, on="case_id", validate="one_to_one")
    donor_table.to_csv(out / "h1_donor_features.csv", index=False)
    for audit in audits:
        source_hashes.update(audit["source_sha256"])
        source_hashes.update(audit["drilldown_sha256"])
        expected = census.loc[census.case_id.eq(audit["case_id"])].iloc[0]
        if audit["analysed_cells_all"] != int(expected.analysed):
            raise ValueError("Unique tissue cell census mismatch")
        if not np.isclose(audit["tissue_area_mm2"], expected.tissue_area_mm2, rtol=1e-9):
            raise ValueError("Tissue area mismatch")
    schema = [{"column": c, "dtype": str(result[c].dtype), "missing": int(result[c].isna().sum())}
              for c in result]
    pd.DataFrame(schema).to_csv(out / "feature_schema.csv", index=False)
    write_json(out / "feature_audit.json", dict(data_root=str(root), full_cohort=not bool(args.cases),
              eligibility_root=str(args.eligibility_root.resolve()) if args.eligibility_root else None,
              cases=cases, objects=len(result), donors=len(cases), donor_groups=groups,
              total_primary_core_cells=int(result.n_cells_core.sum()),
              total_primary_peri_memberships=int(result.n_cells_peri.sum()),
              count_total_scope='all paired source objects, including withheld objects; filter core_complete and neighborhood support for analysis',
              zero_immune_pool_objects=int((~result.immune_exposed).sum()),
              zero_exclusive_immune_pool_objects=int((~result.immune_exposed_exclusive).sum()),
              core_geometry_identity_discrepancies=sum(len(x["core_geometry_identity_discrepancies"]) for x in audits),
              count_discrepancies=sum(len(x["count_discrepancies"]) for x in audits),
              drilldowns_verified=not args.skip_drilldown_verification,
              definitions={"immune_union": sorted(IMMUNE), "endocrine_union": sorted(ENDOCRINE),
                 "primary_peri": "existing tissue-clipped 20um band subtracting every Islet/Cluster; overlap retained",
                 "exclusive_peri": "each eligible actual peri-band cell assigned to nearest eligible Islet/Cluster boundary; stable object-name tie break; no Singles",
                 "local_fields": "20-100um and100-250um boundary annuli clipped to Tissue, subtracting every Islet/Cluster core+20um band; unique cells within each region",
                 "density_units": "cells/mm2 of geometrically assessable region",
                 "zero_pool": "immune localization fraction is missing when core+peri immune count is zero",
                 "cell_identity": "case_id+object_id; never row-order joins",
                 "tissue_background": "unique analysed donor cells with centroids inside saved Tissue divided by donor Tissue area; all analysed counts and outside-Tissue counts retained separately; distal companion excludes all core+20um footprints"},
              software={"python": platform.python_version(), "anndata": version("anndata"), "shapely": shapely.__version__,
                        "numpy": np.__version__, "pandas": pd.__version__},
              script_sha256=sha256(Path(__file__)), argv=sys.argv, source_sha256=source_hashes,
              output_sha256={str(out / f): sha256(out / f) for f in
                             ["h1_object_features.csv", "h1_object_features.parquet", "h1_donor_features.csv", "feature_schema.csv"]},
              donor_audit_files=[str(out / "donor_audits" / f"{case}.json") for case in cases]))
    print(f"Complete: {len(result)} objects/{len(cases)} donors. {out / 'h1_object_features.csv'}", flush=True)


if __name__ == "__main__":
    main()
