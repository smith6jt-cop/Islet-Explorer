#!/usr/bin/env python3
"""Validate and archive final figure/data artifacts without rerunning analyses."""
from pathlib import Path
import hashlib,json,shutil,zipfile
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/app_hypotheses_20261026'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    shutil.copy2(__file__,OUT/'source'/Path(__file__).name)
    figure_index=pd.read_csv(OUT/'figure_index.csv')
    assert len(figure_index)==51
    for column in ['png','pdf','svg','legend']:
        assert all((OUT/p).is_file() for p in figure_index[column])
    sources=json.loads((OUT/'prism/source_hashes.json').read_text())
    assert all(sha(OUT/p)==digest for p,digest in sources.items())
    prism=json.loads((OUT/'prism/validation.json').read_text())
    assert prism['numeric_roundtrip_exact'] and prism['row_identity_roundtrip_exact']
    assert prism['figures']==51 and prism['all_ten_hypotheses_present']
    assert sha(OUT/'prism/App_Hypotheses_plot_data.pzfx')==prism['prism_sha256']
    pdf=json.loads((OUT/'pdf_layout_validation.json').read_text())
    assert pdf['all_51_figures_one_per_page'] and not pdf['text_outside_page']
    assert sha(OUT/'App_Hypotheses_Figures_Statistics.pdf')==pdf['pdf_sha256']
    tests=pd.read_csv(OUT/'primary_statistics.csv')
    assert tests.hypothesis.nunique()==10 and tests.p_value.notna().sum()==9
    execution=pd.read_csv(OUT/'original_plan_execution_audit.csv')
    assert execution.hypothesis.nunique()==10
    qa={'figures':51,'pdf_pages':pdf['pages'],'primary_comparisons':10,
        'available_primary_p_values':9,'prism_native_tables':prism['tables'],
        'prism_numeric_values':prism['numeric_values'],
        'prism_sources_current':len(sources),'conditional_execution_records':len(execution),
        'figure_visual_review':'All 51 figure sets inspected; final 69-page PDF rendered and inspected.',
        'numeric_and_identity_roundtrips':True,'quantitative_INS_endpoints':0,
        'app_pipeline_rebuilt':False,'biological_interpretations_added':False}
    (OUT/'package_validation.json').write_text(json.dumps(qa,indent=2)+'\n')
    files=sorted(p for p in OUT.rglob('*') if p.is_file()
                 and not any(part.startswith('.') or part=='__pycache__' for part in p.relative_to(OUT).parts)
                 and p.name not in {'package_manifest.json','delivery_validation.json'})
    manifest={str(p.relative_to(OUT)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in files}
    (OUT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive=ROOT/'output/App_Hypotheses_20donors_20261026.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[OUT/'package_manifest.json']:
            z.write(p,Path('App_Hypotheses_20donors')/p.relative_to(OUT))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(files)+1
    delivery={'archive':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),
              'files':len(files)+1,'zip_crc_validation':True,'package_checks':qa}
    (OUT/'delivery_validation.json').write_text(json.dumps(delivery,indent=2)+'\n')
    print(json.dumps(delivery,indent=2))

if __name__=='__main__':main()
