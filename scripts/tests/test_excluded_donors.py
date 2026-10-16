"""Donor 78 (HDL078 / HPAP-078) stays out of the app and every cohort loader."""
from pathlib import Path
import re
import sys

import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
REPO = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'research'))
sys.path.insert(0, str(SCRIPTS / 'cohort_expansion'))

from excluded_donors import assert_no_excluded, canonical_id, excluded_ids, find_excluded


@pytest.mark.parametrize('value', ['78', 78, 78.0, '78.0', '0078', 'HPAP-078', 'HPAP078', 'HDL078', 'hdl78', ' 078 '])
def test_every_spelling_of_donor_78_is_excluded(value):
    assert canonical_id(value) == '78'
    with pytest.raises(ValueError, match='excluded donor'):
        assert_no_excluded(['6232', value], 'Cohort')


@pytest.mark.parametrize('value', ['112', 'HDL112', '6232', '780', '1078', '7', 'Islet_78'])
def test_other_donors_pass(value):
    assert find_excluded([value]) == []


def test_registry_lists_only_donor_78():
    assert excluded_ids() == {'78'}


def test_current_protein_additions_and_manifest_exclude_78():
    from stage_original_gate_calls import DONORS
    assert DONORS == ('6090', '6232', '6520', '6539', '6558')
    cohort = pd.read_csv(SCRIPTS / 'protein_expansion_cohort.csv', dtype=str)
    assert find_excluded(cohort.case_id) == []
    assert len(cohort) == 20 and cohort.case_id.is_unique
    assert cohort.donor_status.value_counts().to_dict() == {'ND': 7, 'Aab+': 7, 'T1D': 6}
    assert set(DONORS) <= set(cohort.case_id)


def test_app_constant_matches_registry():
    source = (REPO / 'app/shiny_app/R/00_globals.R').read_text()
    match = re.search(r'^ISLET_EXCLUDED_DONORS <- c\(([^)]*)\)', source, re.M)
    assert match, 'ISLET_EXCLUDED_DONORS is missing from 00_globals.R'
    assert set(re.findall(r'"([^"]+)"', match.group(1))) == excluded_ids()


def test_historical_21_donor_manifest_is_refused(tmp_path):
    from common import read_cohort
    with pytest.raises(ValueError, match='excluded donor'):
        read_cohort(SCRIPTS / 'research/cohort-v1.csv', REPO)


def test_cohort_expansion_manifest_marks_78_excluded(tmp_path):
    import prepare_analysis
    rows = {row['case_id']: row for row in prepare_analysis.read_cohort(SCRIPTS / 'cohort_expansion/cohort.csv')}
    assert rows['78']['inclusion_status'] == 'excluded'
    unmarked = tmp_path / 'cohort.csv'
    unmarked.write_text('case_id,donor_source,source_filename,qupath_entry_id,inclusion_status\n'
                        '78,HPAP,hp.tif,42,selected\n')
    with pytest.raises(ValueError, match='must be marked excluded'):
        prepare_analysis.read_cohort(unmarked)


def test_completed_runners_no_longer_offer_78():
    for name in ['run_qupath.py', 'convert_images.py', 'prepare_qupath.py', 'resume_gpu.py']:
        text = (SCRIPTS / 'cohort_expansion' / name).read_text()
        assert '0078' not in text and 'HPAP-078' not in text, name
