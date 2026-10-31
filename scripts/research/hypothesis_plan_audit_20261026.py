#!/usr/bin/env python3
"""Record plan execution and numerical test inventory without interpreting results."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/h1_h10_conditional_audit_20261026'
REST = ROOT / 'output/h5_h10_original_plan_20261026'
NEW = ROOT / 'output/h1_h4_disease_comparisons_20261026'
PLAN = ROOT / 'docs/T1D_research_hypotheses_20261022.md'

# Execution descriptions, not assessments of the biological hypotheses.
EXECUTION = {
 'H1': ('Count/zero components, corrected immune union and edge exclusions; current disease comparison uses adjusted donor coefficients. Supporting within-state coefficients are reported separately.',
        'The current primary presentation follows the later request for disease comparisons; it is not the original single conditional-count coefficient.',
        'Count, presence/positive-abundance, reciprocal localization, shared/exclusive peri accounting, a 100–250 µm field sensitivity and T1D donor omission checks exist.',
        'The complete 0–20, 20–100 and 100–250 µm comparison and explicit discordant-object inventory are not established as complete by this audit.'),
 'H2': ('Nonvascular ColIV boundary coverage and adjusted core immune density measured from existing images; donor coefficients and disease comparisons exported.',
        'Density regression/donor-coefficient comparisons implement the later conventional-method direction; original count-offset specification is not the current primary model.',
        'Boundary-width sensitivity and computational structural examples exist.',
        'Gap-centered immune distances versus randomized positions on the same perimeter and the complete raw/corrected annulus comparison are not established as complete.'),
 'H3': ('Donor near-crossing minus reference immune-density measurements, paired tests and disease comparisons are exported.',
        'The current comparison uses measured reference compartments; it is not the originally specified geometry-preserving spatial-null statistic.',
        'Crossing-width 5/10/20 µm measurements and comparisons exist.',
        'The current analysis does not complete the separate PDPN-supported lymphatic model or all original structural follow-ups.'),
 'H4': ('Conditional LGALS3 source-allocation model and disease-by-source comparison use accepted cell sources and existing signal measurements.',
        'Later computational source checks supersede earlier proposed manual-review prerequisites; disease contrasts are shown with conditional estimates.',
        'Source thresholds, source-cell overlap, immune localization and alternative eligible populations were checked.',
        'Full composition-matched core/peri/distal comparisons and every secondary source marker are not established as complete.'),
 'H5': ('Additive INS identity, typed-Beta mean INS association, separate hormone-negative/residual contributions and disease comparisons are exported.',
        'Disease comparison is the primary presentation; supporting state-specific conditional means are separately identified.',
        'Raw/corrected INS, median endpoint, strong hormone margin and residual/CHGA accounting are available.',
        'A complete hormone-threshold sweep and an independently defined endocrine-mask analysis are not established as complete.'),
 'H6': ('Combined Beta–Alpha/Delta adjacency minus the label-permutation expectation is tested with composition, area and immune covariates. Original main-association and disease-comparison outputs are separate.',
        'Approved original label permutations are retained. Random-intercept sensitivity is separate from the donor-slope/random-slope specification. Supplementary excess-adjacency group tests are outside the original conditional endpoint.',
        'True polygon-boundary contacts, separate Alpha/Delta, Beta-density endpoint, contact-by-immune interaction, composition/size strata, raw/corrected INS and nerve/vessel adjustment are available.',
        'Independent fully adjusted donor coefficients are estimable for 15 donors, including one T1D donor. The saved primary donor disease test has no p value; the random-intercept sensitivity has a separate result.'),
 'H7': ('Observed INS-independent B3TUBB–CD31 proximity is modeled against within-Beta INS with separate nerve/vessel abundance, size, perimeter and immune covariates.',
        'The user-approved observed-proximity regression replaces excess-over-null alignment. No excess-over-null result is asserted. Random-intercept sensitivity is separately identified.',
        '5/10/20 µm proximity, overlap, GAP43, raw/corrected INS and abundance/geometry companion measurements are available.',
        'Core/distal fiber-area measurements, fiber-length comparisons and all inherited-mask comparisons were not run; peri area was measured. Fully adjusted primary individual coefficients are estimable for 16 donors, including two T1D donors.'),
 'H8': ('Full geometric census including unpaired objects and object-free tissue tiles; 2–9 typed endocrine-cell Cluster density modeled against fields, with disease comparisons.',
        'Conventional donor-coefficient comparisons implement the disease-by-field analysis; census and eligibility definitions are explicitly retained.',
        '100/200/400 µm grids, raw whole-cell count, geometric size, uncertain-identity bounds and object-class density companions are available.',
        'Age-matched subsets and the optional Ki67 readout were not run.'),
 'H9': ('Original shared-range requirement checked: one Aab+ and two T1D donors occupy the common Beta-quantity range; a three-parameter model there has zero residual degrees of freedom.',
        'Original primary p is unavailable. The 13-donor all-range ANCOVA remains a separately labeled supplementary calculation; it does not meet the original shared-range restriction.',
        'Shared-range donor and object counts plus age and regional-immune adjusted supplementary models are exported.',
        'Single/multiple-Aab and antibody-specific neighborhood follow-ups are not established as complete.'),
 'H10': ('Approved nested log-C-peptide models and donor-held-out RMSE use 19 measured donors and the fixed score with training-only scaling. Proper censored likelihood uses all 20 donors, including the <0.02 observation.',
         'Ordinary nested regression and RMSE replace regularized predictive log scores/conditional simulation. Full-sample F/LR tests and held-out error are separate outputs.',
         'Disease-conditioned models, fixed-score H5 INS prediction, disease-specific errors and all six individual T1D observations/predictions are exported.',
         'Additional HbA1c, peak glucose, BMI and procurement sensitivity models were not run. No T1D-only predictive model is substituted for the approved analysis.'),
}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    text = PLAN.read_text()
    sections = re.split(r'(?m)^### (H\d+) ', text)
    original = {}
    for i in range(1, len(sections), 2):
        h, body = sections[i], sections[i + 1]
        fields = {}
        for label, content in re.findall(r'\*\*([^*]+)\.\*\*\s*(.*?)(?=\n\n|\Z)', body, flags=re.S):
            if label in ['Hypothesis', 'Primary test', 'If supported', 'Same-data follow-up'] or label.startswith('If a meaningful'):
                fields[label] = content.strip()
        original[h] = fields
    assert len(original) == 10 and all('If supported' in x for x in original.values())
    records = []
    for h, (execution, adaptation, completed, remaining) in EXECUTION.items():
        negative = next(k for k in original[h] if k.startswith('If a meaningful'))
        records.append(dict(hypothesis=h, execution=execution, adaptation_or_difference=adaptation,
            followups_available=completed, remaining_or_technical_limits=remaining,
            original_if_supported=original[h]['If supported'], original_exclusion_condition=negative,
            original_if_excluded=original[h][negative], result_branch='Not assigned',
            reason='Owner requested no autonomous result interpretation; this is an execution/compliance audit.'))
    pd.DataFrame(records).to_csv(OUT / 'plan_execution_audit.csv', index=False)
    (OUT / 'original_plan_conditions.json').write_text(json.dumps(original, indent=2, ensure_ascii=False) + '\n')

    tests = pd.read_csv(NEW / 'primary_disease_comparison_tests.csv')[['hypothesis', 'relationship', 'test', 'p']]
    more = [
        ('H5', 'Disease comparison of immune–INS coefficients', 'Kruskal–Wallis', pd.read_csv(REST/'h5/donor_disease_global_test.csv').p_value.iloc[0], 'Computed'),
        ('H6', 'Disease comparison of adjusted adjacency–INS coefficients', 'Independent donor-coefficient comparison', np.nan, 'No primary p: one estimable T1D donor coefficient; separate random-intercept sensitivity retained'),
        ('H7', 'Disease comparison of observed proximity–INS coefficients', 'Kruskal–Wallis', pd.read_csv(REST/'h7/fallback_disease_tests.csv').query("model == 'primary_20um'").p.iloc[0], 'Computed'),
        ('H8', 'Disease comparison of immune-field small-Cluster coefficients', 'Kruskal–Wallis', pd.read_csv(REST/'h8/disease_field_overall_tests.csv').query("model == 'primary_200um' and field == 'immune'").p.iloc[0], 'Computed'),
        ('H9', 'Shared-range Aab+–T1D coverage association', 'Donor ANCOVA within shared Beta range', np.nan, 'Zero residual degrees of freedom in shared range'),
        ('H10', 'Added fixed-score association with log C-peptide', 'Nested partial F test', pd.read_csv(REST/'h10/primary_test.csv').p_value.iloc[0], 'Computed; held-out RMSE reported separately'),
    ]
    tests['execution_status'] = 'Computed'
    tests = pd.concat([tests, pd.DataFrame(more, columns=tests.columns)], ignore_index=True)
    valid = tests.p.notna()
    ranked = tests.loc[valid, 'p'].sort_values()
    # BH over the original ten declared slots; unavailable endpoints stay missing.
    adjusted = np.minimum.accumulate((ranked.to_numpy() * 10 / np.arange(1, len(ranked)+1))[::-1])[::-1]
    tests['p_BH_ten_slots_current_comparisons'] = np.nan
    tests.loc[ranked.index, 'p_BH_ten_slots_current_comparisons'] = np.minimum(1, adjusted)
    tests['family_scope'] = 'Current disease comparisons H1–H8, shared-range H9, approved nested H10; not the original ten conditional-direction tests'
    tests.to_csv(OUT / 'current_test_inventory.csv', index=False)
    assert len(tests) == 10 and tests.loc[tests.hypothesis.isin(['H6','H9']), 'p'].isna().all()

    lines = ['# Original hypothesis-plan execution audit', '',
        'This document records methods, completed work, numerical test locations and departures from the original plan. '
        'It does not interpret the results or assign biological support/exclusion branches. Earlier agent-authored conclusions are not adopted by this audit.', '',
        'The original plan is preserved at `docs/T1D_research_hypotheses_20261022.md`. Later owner instructions govern conventional methods, '
        'computational source validation, primary disease comparisons, and the approved H6/H7/H10 methods. The latest instruction restricts interpretation.', '',
        '## Conditional-result rules', '',
        '- Disease-state contrasts and conditional exposure–outcome associations are separately labeled; they test different coefficients.',
        '- The original negative branches require exclusion of a meaningful effect, not p above 0.05. No new meaningful-effect threshold is introduced here.',
        '- Original “If supported” and “If a meaningful effect is excluded” text is retained in `original_plan_conditions.json`, without assigning either branch.',
        '- H9 explicitly specifies an unestimable status when overlap/replication requirements fail. Its shared-range model has zero residual degrees of freedom. The all-range model is supplementary.',
        '- H10 association tests and prediction errors are kept separate. The approved change in method is recorded rather than described as an exact execution of the original predictive-log-score specification.',
        '- Conditional same-data follow-ups are inventoried below. The complete list of optional follow-ups has not been executed.', '',
        '## Execution and differences', '',
        '| Hypothesis | Executed | Adaptation or difference | Available follow-ups | Remaining work or technical limits |',
        '|---|---|---|---|---|']
    for r in records:
        lines.append('| ' + ' | '.join(str(r[k]).replace('|','/') for k in ['hypothesis','execution','adaptation_or_difference','followups_available','remaining_or_technical_limits']) + ' |')
    lines += ['', '## Numerical test inventory', '',
        'The table follows the current requested disease-comparison presentation and approved H10 nested test. '
        'It is not a table of the original ten conditional-direction tests. BH adjustment reserves all ten slots, '
        'retains unavailable p values as missing and does not substitute supporting sensitivity models for unavailable primary tests. '
        'Each figure’s previously declared pairwise or sensitivity family remains separately documented in its legend.', '',
        '| Hypothesis | Current comparison p | BH across ten slots | Execution status |', '|---|---:|---:|---|']
    for r in tests.itertuples():
        p = 'Unavailable' if pd.isna(r.p) else f'{r.p:.6g}'
        q = 'Unavailable' if pd.isna(r.p_BH_ten_slots_current_comparisons) else f'{r.p_BH_ten_slots_current_comparisons:.6g}'
        lines.append(f'| {r.hypothesis} | {p} | {q} | {r.execution_status} |')
    lines += ['', 'Supporting H1/H2/H5 state-specific estimates, intervals and numerical tests are in `association_audit/`. '
        'H6–H10 results, specifications and exact plot data are in `output/h5_h10_original_plan_20261026/`. '
        'The shared Prism package includes these results and the source test tables. No new biological interpretation or follow-up recommendation is made here.', '']
    (OUT / 'README.md').write_text('\n'.join(lines))
    provenance = dict(plan=str(PLAN.relative_to(ROOT)), plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        hypotheses=10, branch_assignments=0, numerical_primary_p_values=int(valid.sum()), declared_BH_slots=10,
        missing_primary_p_hypotheses=tests.loc[~valid,'hypothesis'].tolist(), original_plan_modified=False,
        owner_interpretation_restriction_applied=True)
    (OUT / 'validation.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance))

if __name__ == '__main__':
    main()
