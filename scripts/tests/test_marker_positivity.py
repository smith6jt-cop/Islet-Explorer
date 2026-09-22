"""Nullable, evidence-bound marker calls preserve excluded raw objects."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from common import MARKERS,sha256,write_json
from marker_positivity import apply,reviewed_specification


def test_unreviewed_markers_and_segmentation_fragments_never_become_negative(tmp_path):
    root=tmp_path/'candidate';out=root/'positivity';out.mkdir(parents=True)
    corrected=root/'corrected';corrected.mkdir()
    frame=pd.DataFrame(dict(case_id=['test']*3,object_id=['test:1','test:2','test:3'],
                            **{'Cell Area':[20,20,.2]},x_source_px=[0,20,30],y_source_px=[0,20,30]))
    for m in MARKERS:frame[m]=[10,100,1000]
    path=corrected/'test-cells.parquet';frame.to_parquet(path,index=False)
    write_json(out/'candidate_specification.json',dict(input_hashes={'test':sha256(path)},technical_exclusion_regions=[]))
    evidence=tmp_path/'synthetic_review.txt';evidence.write_text('Synthetic test evidence, not biological validation.')
    records={m:dict(validated=False,threshold=None,evidence=[],reason='Unreviewed') for m in MARKERS}
    records['INS']=dict(validated=True,threshold=50,reason='Synthetic gate',evidence=[dict(path=str(evidence),sha256=sha256(evidence))])
    review=out/'review.json';write_json(review,dict(reviewer='Unit fixture',markers=records,candidate_specification_sha256=sha256(out/'candidate_specification.json')))
    apply(root,out,review);calls=pd.read_parquet(out/'test.parquet')
    assert calls.object_id.tolist()==frame.object_id.tolist()
    assert calls.positive__INS.iloc[:2].tolist()==[False,True]
    assert pd.isna(calls.positive__INS.iloc[2]) and calls.positive__CD45.isna().all()
    evidence.write_text('Changed evidence')
    with pytest.raises(ValueError,match='evidence changed'):reviewed_specification(out,review)
