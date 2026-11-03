#!/usr/bin/env python3
"""Compose/validate the narrowly scoped, evidence-bound cell-state revalidation.

The failed full report remains immutable. Its runtime/producer identity is never
replaced with the separately versioned corrected validator identity.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

PRODUCER = 'e5689f6453e68bc9e94cbe83689cb9db4d3d0b83'
VALIDATOR = 'd3ec00e113fbbdb8abb7eb4b57ee4ffc7ca46580'
VERIFIER_SHA = '672658b045514616b1b32bd65afa4be5904d7a3778a5951157558e1b566dc648'
SCHEMA = 'islet.parity_targeted_revalidation/1'
METHOD = '24 unchanged full-run checks plus corrected complete cell-state gate and independent UUID audit'
CHECKS = {'canonical_h5ad_preparation', 'release_manifest_and_full_hashes', 'release_source_identity',
          'canonical_source_hashes', 'endocrine_group_definition', 'passthrough_files_and_trees',
          'embedding_metadata', 'signed_distance_kde_table', 'donor_tissue_tables', 'cell_state_exports',
          'spatial_source_load', 'geometry_representative_bounding_boxes'}
CHECKS.update('summary_' + x for x in ('core_area', 'targets_all', 'markers_all', 'comp', 'combined_pairs'))
CHECKS.update(f'spatial_{kind}_{mode}_{norm}' for kind in ('histogram', 'composition')
              for mode in ('phenotype', 'marker') for norm in ('prop', 'count'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def reference(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': sha(path)}


def load_reference(ref):
    require(set(ref) == {'path', 'sha256'}, 'invalid evidence reference')
    require(sha(ref['path']) == ref['sha256'], 'evidence hash mismatch')
    return read(ref['path'])


def identity(path):
    path = Path(path)
    st = path.stat()
    return {'path': str(path.resolve(strict=True)), 'device': st.st_dev, 'inode': st.st_ino,
            'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns, 'ctime_ns': st.st_ctime_ns}


def link_identity(path):
    st = path.lstat()
    return {'target': os.readlink(path), 'device': st.st_dev, 'inode': st.st_ino,
            'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns, 'ctime_ns': st.st_ctime_ns}


def relative_path(root, relative):
    p = Path(relative)
    require(not p.is_absolute() and '..' not in p.parts, 'unsafe relative path')
    result = root / p
    require(result.resolve().is_relative_to(root.resolve()), 'path escaped root')
    return result


def tree_files(root):
    return sorted((p for p in root.rglob('*') if p.is_file()), key=lambda p: p.relative_to(root).as_posix().encode())


def digest_tree(records):
    return hashlib.sha256('\n'.join(k + '  ' + v for k, v in sorted(records.items(), key=lambda x: x[0].encode())).encode()).hexdigest()


def check_original(original):
    require(original['schema_version'] == 1 and original['mode'] == 'release_comparison'
            and original['status'] == 'failed', 'not the failed full comparison')
    require(original['numeric_absolute_tolerance'] == 1e-12
            and original['categorical_count_na_parity'] == 'exact', 'comparison semantics changed')
    ident = original['identity']
    require(ident['git_commit'] == ident['release_git_commit'] == PRODUCER and ident['git_dirty'] is False,
            'original producer identity mismatch')
    require(len(original['checks']) == 25 and {r['check'] for r in original['checks']} == CHECKS,
            'incomplete original check inventory')
    failures = [r for r in original['checks'] if r['status'] != 'passed']
    require(failures == [{'check': 'cell_state_exports', 'status': 'failed',
                         'details': {'error_code': 'numeric_mismatch'}}], 'original sole failure mismatch')


def check_correction(correction, uuid, original, source_manifest, release_manifest):
    require(correction['schema'] == 'uniform.verifier_correction_realdata/1' and correction['passed'] is True,
            'correction did not pass')
    require(correction['verifier_sha256'] == VERIFIER_SHA, 'wrong corrected verifier')
    require(correction['release_manifest_sha256'] == original['identity']['release_manifest_sha256'], 'release binding mismatch')
    require(correction['negative_controls'] == {'binary_flag_change': 'numeric_mismatch',
            'distance_change_over_tolerance': 'numeric_mismatch'}, 'missing corruption rejection controls')
    details = correction['complete_cell_state_exports']
    require(details == {'applicable': True, 'fields': ['Ki67_pos'], 'distance_rows': 3199606,
            'donor_files': 20, 'donor_rows': 3199606, 'object_files': 4848, 'object_rows': 357658}, 'incomplete corrected state gate')
    distance = correction['distance_comparison']
    require(set(distance) == {'rows', 'columns', 'keys', 'max_absolute_difference'}
            and distance['rows'] == 3199606 and distance['columns'] == 6
            and distance['keys'] == ['case_id', 'donor_status', 'cell_region', 'phenotype'],
            'incomplete corrected distance comparison')
    require(math.isfinite(distance['max_absolute_difference']) and 0 <= distance['max_absolute_difference'] <= 1e-12,
            'distance exceeds existing tolerance')
    require(uuid['schema'] == 'uniform.independent_uuid_cell_state_release_check/1' and uuid['result'] == 'pass', 'UUID audit failed')
    require(uuid['cells'] == 3199606 and uuid['donors'] == 20 and uuid['ki67_positive_cells'] == 78680
            and uuid['ki67_mismatches'] == 0 and uuid['source_and_release_distance_row_order_exact_for_all_categorical_fields_and_state'] is True,
            'UUID/state mismatch')
    rows = uuid['donor_rows']
    require(len(rows) == len({r['case_id'] for r in rows}) == 20 and sum(r['cells'] for r in rows) == 3199606,
            'UUID donor inventory mismatch')
    require(all(r['cells'] > 0 and r['cells'] == r['unique_uuid_count'] and r['ki67_mismatches'] == 0
            and r['all_source_and_release_uuid_row_associations_exact'] is True
            and set(r['coordinate_comparison']) == {'X_centroid', 'Y_centroid'}
            and all(c['differing'] == 0 and c['max_absolute_difference'] == 0 for c in r['coordinate_comparison'].values()) for r in rows), 'UUID row association mismatch')
    require(sum(r['ki67_positive_cells'] for r in rows) == 78680, 'UUID positive-state total mismatch')
    ud = uuid['distance_measurement']
    require(set(ud) == {'csv_reader', 'parquet_reader', 'differing', 'max_absolute_difference', 'differences_over_existing_1e_12_tolerance'}
            and isinstance(ud['differing'], int) and 0 <= ud['differing'] <= 3199606
            and math.isfinite(ud['max_absolute_difference']) and 0 <= ud['max_absolute_difference'] <= 1e-12
            and ud['differences_over_existing_1e_12_tolerance'] == 0, 'incomplete UUID distance evidence')
    require(uuid['object_csv_files_byte_identical'] == 4848, 'object passthrough mismatch')
    require(uuid['source_distance_sha256'] == correction['source_cell_distance_sha256'], 'distance source binding mismatch')
    require(uuid['source'] == correction['source_root'] and uuid['release'] == correction['release_root'], 'audit path mismatch')
    require(release_manifest['git_commit'] == PRODUCER and release_manifest['git_dirty'] is False
            and str(release_manifest['schema_version']) == '5', 'release producer changed')
    require(source_manifest['git_commit'] == '67ee6c9d0b28f467eee589fea24829a6546e3d74', 'scientific producer changed')


def capture_inputs(source, release, original, inventory, image_expectation):
    """Rehash small inputs; retain original full image proof only before its start."""
    cutoff = int(datetime.fromisoformat(original['created_at_utc'].replace('Z', '+00:00')).timestamp() * 1e9)
    release_manifest = read(release / 'manifest.json')
    records, rosters, links = {}, {}, {}

    def capture(path, expected, hash_now=True):
        before = identity(path)
        require(before['ctime_ns'] <= cutoff and before['mtime_ns'] <= cutoff,
                'input modified after original full verification began')
        if hash_now:
            require(sha(path) == expected, 'input content hash changed')
        require(identity(path) == before, 'input changed while checking')
        records[str(path)] = {'identity': before, 'sha256': expected,
                              'method': 'fresh_sha256' if hash_now else 'original_full_hash_plus_unchanged_inode_times'}

    source_files = inventory['roots']['app']['files']
    for relative, row in source_files.items():
        p = relative_path(source, relative)
        require(p.stat().st_size == row['bytes'], 'source input size mismatch')
        capture(p, row['sha256'])
    actual_source = {p.relative_to(source).as_posix() for p in tree_files(source) if not p.is_symlink()}
    require(actual_source == set(source_files), 'source input roster changed')
    rosters[str(source)] = sorted(source_files)
    for row in release_manifest['canonical_sources']:
        require(source_files[row['path']]['sha256'] == row['sha256'], 'canonical source hash mismatch')
    image_rows = {r['relative_path']: r for r in image_expectation['files']}
    for relative, target in inventory['roots']['app']['symlinks'].items():
        path = source / relative
        require(path.is_symlink() and os.readlink(path) == target, 'source image link changed')
        link, resolved = link_identity(path), identity(path)
        require(all(s['mtime_ns'] <= cutoff and s['ctime_ns'] <= cutoff for s in (link, resolved)), 'source image changed after full verification began')
        expected = image_rows[Path(relative).name]
        require(resolved['bytes'] == expected['bytes'], 'source image size changed')
        links[relative] = {'link': link, 'target': resolved, 'sha256': expected['sha256']}
    require(image_expectation['sha256'] == release_manifest['artifacts']['images']['sha256'], 'image expectation changed')
    for name, artifact in release_manifest['artifacts'].items():
        path = relative_path(release, artifact['relative_path'])
        if artifact['type'] == 'file':
            capture(path, artifact['sha256'])
        else:
            files = tree_files(path)
            require(len(files) == artifact['file_count'], 'release tree inventory changed')
            hashes = {}
            for p in files:
                rel = p.relative_to(path).as_posix()
                require(not p.is_symlink(), 'unexpected release symlink')
                if name == 'images':
                    require(rel in image_rows and p.stat().st_size == image_rows[rel]['bytes'], 'image roster or size changed')
                    digest = image_rows[rel]['sha256']
                    capture(p, digest, hash_now=False)
                else:
                    digest = sha(p)
                    capture(p, digest)
                hashes[rel] = digest
            require(digest_tree(hashes) == artifact['sha256'], 'release tree hash changed')
            rosters[str(path)] = sorted(hashes)
    require(sha(source / 'islet_explorer.h5ad') == original['identity']['canonical_h5ad_sha256'], 'H5AD changed')
    require(sha(source / 'donor_metadata.csv') == original['identity']['donor_metadata_sha256'], 'donor metadata changed')
    return {'schema': 'islet.parity_unchanged_inputs/1', 'original_started_utc': original['created_at_utc'],
            'source_root': str(source), 'release_root': str(release), 'files': records, 'rosters': rosters, 'source_links': links,
            'source_inventory': reference(inventory_path_global), 'image_expectation': reference(image_expectation_path_global),
            'image_bytes_rehashed': 0, 'scope': 'Original full hashes reused only for the exact image roster with unchanged identity and pre-verification ctime/mtime; all non-image source and release content hashes rechecked.'}


def validate_inputs(proof, manifest):
    require(proof['schema'] == 'islet.parity_unchanged_inputs/1' and proof['image_bytes_rehashed'] == 0, 'input proof schema mismatch')
    inventory = load_reference(proof['source_inventory'])
    images = load_reference(proof['image_expectation'])
    source, release = Path(proof['source_root']), Path(proof['release_root'])
    cutoff = int(datetime.fromisoformat(proof['original_started_utc'].replace('Z', '+00:00')).timestamp() * 1e9)
    expected, rosters = {}, {str(source): sorted(inventory['roots']['app']['files'])}
    for relative, row in inventory['roots']['app']['files'].items():
        expected[str(relative_path(source, relative))] = (row['sha256'], row['bytes'], 'fresh_sha256')
    image_rows = {r['relative_path']: r for r in images['files']}
    require(images['source_root'] == str(source) and images['sha256'] == manifest['artifacts']['images']['sha256'], 'image evidence root/hash mismatch')
    for name, artifact in manifest['artifacts'].items():
        path = relative_path(release, artifact['relative_path'])
        if artifact['type'] == 'file':
            expected[str(path)] = (artifact['sha256'], artifact['bytes'], 'fresh_sha256')
            continue
        files = tree_files(path)
        require(len(files) == artifact['file_count'] and sum(p.stat().st_size for p in files) == artifact['bytes'], 'release tree count/size mismatch')
        hashes = {}
        for p in files:
            require(not p.is_symlink(), 'release tree symlink')
            rel = p.relative_to(path).as_posix()
            require(str(p) in proof['files'], 'omitted input proof entry')
            if name == 'images':
                require(rel in image_rows, 'image roster mismatch')
                digest, size = image_rows[rel]['sha256'], image_rows[rel]['bytes']
                method = 'original_full_hash_plus_unchanged_inode_times'
            else:
                digest, size, method = proof['files'][str(p)]['sha256'], p.stat().st_size, 'fresh_sha256'
            expected[str(p)] = (digest, size, method)
            hashes[rel] = digest
        require(digest_tree(hashes) == artifact['sha256'], 'proof tree hash mismatch')
        rosters[str(path)] = sorted(hashes)
    require(set(proof['files']) == set(expected) and proof['rosters'] == rosters, 'incomplete proof entries or rosters')
    expected_release = {str(release / 'manifest.json')} | {p for p in expected if Path(p).is_relative_to(release)}
    require(not any(p.is_symlink() for p in release.rglob('*')) and {str(p) for p in tree_files(release)} == expected_release,
            'unlisted release file or symlink')
    for path, row in proof['files'].items():
        digest, size, method = expected[path]
        require(row['sha256'] == digest and row['method'] == method and row['identity']['bytes'] == size, 'proof hash/method/size mismatch')
        require(row['identity']['ctime_ns'] <= cutoff and row['identity']['mtime_ns'] <= cutoff, 'proof input postdates original verification')
        require(identity(path) == row['identity'], 'input identity changed since revalidation')
    for root, roster in proof['rosters'].items():
        path = Path(root)
        actual = sorted(p.relative_to(path).as_posix() for p in tree_files(path) if not p.is_symlink())
        require(actual == roster, 'input roster changed since revalidation')
    expected_links = inventory['roots']['app']['symlinks']
    require(set(proof['source_links']) == set(expected_links) == {p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_symlink()}, 'source link roster mismatch')
    for relative, target in expected_links.items():
        p, row = source / relative, proof['source_links'][relative]
        require(row['link']['target'] == target and link_identity(p) == row['link'] and identity(p) == row['target'], 'source image link/target changed')
        require(all(s['mtime_ns'] <= cutoff and s['ctime_ns'] <= cutoff for s in (row['link'], row['target'])), 'source link proof postdates full verification')
        require(row['sha256'] == image_rows[Path(relative).name]['sha256'] and row['target']['bytes'] == image_rows[Path(relative).name]['bytes'], 'source image pin mismatch')


def expected_composite(original, correction, provenance, created):
    result = copy.deepcopy(original)
    result['created_at_utc'] = created
    result['status'] = 'passed'
    result['checks'] = [({'check': 'cell_state_exports', 'status': 'passed',
                         'details': correction['complete_cell_state_exports']} if row['check'] == 'cell_state_exports' else row)
                        for row in original['checks']]
    result['validation_provenance'] = provenance
    return result


def validate_report(path, repo=None):
    report = read(path)
    p = report['validation_provenance']
    require(set(p) == {'schema', 'validator_commit', 'verifier_sha256', 'composer_sha256', 'method', 'evidence'}
            and p['method'] == METHOD and p['composer_sha256'] == sha(Path(__file__)), 'composer/revalidation metadata mismatch')
    require(p['schema'] == SCHEMA and p['validator_commit'] == VALIDATOR and p['verifier_sha256'] == VERIFIER_SHA,
            'validator identity mismatch')
    evidence = {k: load_reference(v) for k, v in p['evidence'].items()}
    require(set(evidence) == {'original_full_report', 'corrected_state_report', 'independent_uuid_audit',
            'source_manifest', 'release_manifest', 'unchanged_inputs'}, 'evidence inventory mismatch')
    original, correction = evidence['original_full_report'], evidence['corrected_state_report']
    check_original(original)
    check_correction(correction, evidence['independent_uuid_audit'], original, evidence['source_manifest'], evidence['release_manifest'])
    require(p['evidence']['release_manifest']['sha256'] == original['identity']['release_manifest_sha256']
            and p['evidence']['source_manifest']['sha256'] == correction['source_uniform_manifest_sha256'], 'manifest binding mismatch')
    proof = evidence['unchanged_inputs']
    require(proof['original_started_utc'] == original['created_at_utc']
            and proof['source_root'] == correction['source_root'] and proof['release_root'] == correction['release_root'], 'input proof identity mismatch')
    require(proof['files'][str(Path(correction['source_root']) / 'cell_distance.csv')]['sha256'] == correction['source_cell_distance_sha256'], 'source state bytes mismatch')
    require(proof['files'][str(Path(correction['release_root']) / 'spatial/cell_distance.parquet')]['sha256'] == evidence['independent_uuid_audit']['release_distance_sha256'], 'release state bytes mismatch')
    validate_inputs(proof, evidence['release_manifest'])
    repo = Path(repo or Path(__file__).resolve().parent.parent)
    blob = subprocess.check_output(['git', '-C', str(repo), 'show', VALIDATOR + ':scripts/verify_release_parity.R'])
    require(hashlib.sha256(blob).hexdigest() == VERIFIER_SHA, 'committed validator blob mismatch')
    require(report == expected_composite(original, correction, p, report['created_at_utc']), 'composite altered retained checks or provenance')
    return {'passed': True, 'report_sha256': sha(path), 'producer_commit': PRODUCER,
            'validator_commit': VALIDATOR, 'reused_full_checks': 24, 'revalidated_checks': ['cell_state_exports']}


def write_new(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('compose')
    for name in ('original', 'correction', 'uuid-audit', 'source-inventory', 'image-expectation', 'input-proof', 'output'):
        build.add_argument('--' + name, type=Path, required=True)
    verify = sub.add_parser('validate'); verify.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'validate':
        print(json.dumps(validate_report(args.report))); return
    require(not args.output.exists() and not args.input_proof.exists(), 'outputs must be fresh')
    original, correction, uuid = read(args.original), read(args.correction), read(args.uuid_audit)
    check_original(original)
    source, release = Path(correction['source_root']), Path(correction['release_root'])
    check_correction(correction, uuid, original, read(source / 'uniform_release.json'), read(release / 'manifest.json'))
    global inventory_path_global, image_expectation_path_global
    inventory_path_global, image_expectation_path_global = args.source_inventory, args.image_expectation
    proof = capture_inputs(source, release, original, read(args.source_inventory), read(args.image_expectation))
    write_new(args.input_proof, proof)
    provenance = {'schema': SCHEMA, 'validator_commit': VALIDATOR, 'verifier_sha256': VERIFIER_SHA,
                  'composer_sha256': sha(Path(__file__)), 'method': METHOD,
                  'evidence': {k: reference(v) for k, v in {
                      'original_full_report': args.original, 'corrected_state_report': args.correction,
                      'independent_uuid_audit': args.uuid_audit, 'source_manifest': source / 'uniform_release.json',
                      'release_manifest': release / 'manifest.json', 'unchanged_inputs': args.input_proof}.items()}}
    composite = expected_composite(original, correction, provenance, datetime.now(timezone.utc).isoformat())
    write_new(args.output, composite)
    print(json.dumps(validate_report(args.output)))


if __name__ == '__main__':
    main()
