#!/usr/bin/env python3
"""Copy the six completed native QuPath hierarchies into an editable review project."""
import argparse
import copy
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlparse


SELECTION = [('78', '', 42), ('6090', 'recovery', 39), ('6232', '', 43),
             ('6520', '', 40), ('6539', 'recovery', 38), ('6558', '', 44)]


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def image_uris(value):
    if isinstance(value, dict):
        for key, child in value.items():
            if key == 'uri':
                yield child
            else:
                yield from image_uris(child)
    elif isinstance(value, list):
        for child in value:
            yield from image_uris(child)


def main():
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path,
                        default=repo/'data/development/expanded-20260911')
    parser.add_argument('--output-root', type=Path,
                        default=repo/'data/reviews/qupath-six-donors-20260917')
    args = parser.parse_args()
    source_root, output = args.source_root.resolve(), args.output_root.resolve()
    if output.exists():
        raise SystemExit(f'Review directory exists; preserve review edits: {output}')
    with (repo/'scripts/research/cohort-v1.csv').open() as stream:
        cohort = {r['case_id']: r for r in csv.DictReader(stream)}
    project = copy.deepcopy(json.loads((source_root/'qupath/project.qpproj').read_text()))
    project['images'] = []
    records = []
    for case_id, suffix, entry_id in SELECTION:
        root = source_root/suffix
        source_project = root/'qupath/project.qpproj'
        document = json.loads(source_project.read_text())
        entry = copy.deepcopy(next(e for e in document['images'] if e['entryID'] == entry_id))
        folder = root/'qupath/data'/str(entry_id)
        summary = json.loads((folder/'summary.json').read_text())
        export_dir = root/'exports'/('0078' if case_id == '78' else case_id)
        exported = json.loads((export_dir/'export_manifest.json').read_text())
        cells = summary['hierarchy']['objectTypeCounts']['Cells']
        if str(exported['case_id']) != case_id or cells != exported['cell_count']:
            raise ValueError(f'Cell export / saved hierarchy mismatch: {case_id}')
        if entry['imageName'] != exported['source_filename']:
            raise ValueError(f'Image name / export mismatch: {case_id}')
        uris = list(image_uris(entry['serverBuilder']))
        if len(uris) != 1 or urlparse(uris[0]).scheme != 'file':
            raise ValueError(f'Expected one local source image: {case_id}')
        source_image = Path(unquote(urlparse(uris[0]).path))
        if not source_image.is_file():
            raise FileNotFoundError(source_image)
        if case_id == '6539' and source_image.resolve().parent != (source_root/'recovery/sources').resolve():
            raise ValueError('6539 must use the recovered acquisition')
        entry.setdefault('metadata', {}).update({
            'Review donor': 'HPAP-078' if case_id == '78' else case_id,
            'Disease status': cohort[case_id]['donor_status'],
            'Annotation review': cohort[case_id]['annotation_review'],
            'Phenotype status': 'Retained subtype assignments await gate validation',
        })
        project['images'].append(entry)
        records.append(dict(case_id=case_id, donor_status=cohort[case_id]['donor_status'],
                            entry_id=entry_id, source_project=str(source_project),
                            source_data=str(folder/'data.qpdata'), source_image=str(source_image),
                            viewer_image=str((repo/cohort[case_id]['image']).resolve()),
                            cell_count=cells, islet_count=summary['hierarchy']['annotationClassificationCounts']['Islet'],
                            export_manifest=str(export_dir/'export_manifest.json'),
                            export_manifest_sha256=digest(export_dir/'export_manifest.json'),
                            annotation_review=cohort[case_id]['annotation_review']))
    output.mkdir(parents=True)
    for record in records:
        original = Path(record['source_data'])
        target = output/'data'/str(record['entry_id'])
        target.mkdir(parents=True)
        for name in ('data.qpdata', 'server.json', 'summary.json'):
            # Independent copies: editing review data must not alter processing evidence.
            shutil.copy2(original.parent/name, target/name)
        copied = target/'data.qpdata'
        if original.stat().st_ino == copied.stat().st_ino and original.stat().st_dev == copied.stat().st_dev:
            raise ValueError('Review objects must not be hardlinked to source objects')
        checksum = digest(original)
        if checksum != digest(copied):
            raise ValueError(f'Copied object data hash mismatch: {record["case_id"]}')
        record.update(review_data=str(copied), data_sha256=checksum, data_bytes=copied.stat().st_size)
        print(f'{record["case_id"]}: copied and verified {record["cell_count"]:,} cells', flush=True)
    classes = source_root/'qupath/classifiers/classes.json'
    (output/'classifiers').mkdir()
    shutil.copy2(classes, output/'classifiers/classes.json')
    project['uri'] = (output/'project.qpproj').as_uri()
    project['modifyTimestamp'] = int(datetime.now(timezone.utc).timestamp()*1000)
    (output/'project.qpproj').write_text(json.dumps(project, indent=2)+'\n')
    manifest = dict(status='copied_and_hash_verified', created_at_utc=datetime.now(timezone.utc).isoformat(),
                    qupath_version=project['version'], project=str(output/'project.qpproj'),
                    project_sha256=digest(output/'project.qpproj'), donor_count=len(records),
                    cell_count=sum(r['cell_count'] for r in records),
                    islet_count=sum(r['islet_count'] for r in records),
                    source_objects_modified=False, image_pixels_copied=False, donors=records)
    (output/'review_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    with (output/'donors.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(output/'project.qpproj', flush=True)


if __name__ == '__main__':
    main()
