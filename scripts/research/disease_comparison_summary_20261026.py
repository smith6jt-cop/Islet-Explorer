#!/usr/bin/env python3
"""Assemble direct disease-state comparisons and their manuscript legends."""
from pathlib import Path
import hashlib
import json
import re

import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/h1_h4_disease_comparisons_20261026'
OLD = ROOT / 'output/h1_h4_literature_methods_20261025'


def main():
    h1 = pd.read_csv(OUT/'h1/donor_disease_global_test.csv').iloc[0]
    h2 = pd.read_csv(OUT/'h2/primary_disease_overall_test.csv').iloc[0]
    h3 = pd.read_csv(OUT/'h3/h3_global_test.csv').iloc[0]
    h4 = pd.read_csv(OUT/'h4/model_anova.csv')
    h4 = h4.loc[h4.model.eq('primary') & h4.term.str.contains(':') & h4.term.str.contains('source_z')].iloc[0]
    summary = pd.DataFrame([
        dict(hypothesis='H1', relationship='Core immune allocation and β-cell density',
             test='Kruskal–Wallis on adjusted donor coefficients', p=float(h1.p_value)),
        dict(hypothesis='H2', relationship='Boundary ColIV coverage and core immune density',
             test='Kruskal–Wallis on adjusted donor coefficients', p=float(h2.p)),
        dict(hypothesis='H3', relationship='Near-capillary minus reference immune density',
             test='Kruskal–Wallis on donor density differences', p=float(h3.p_value)),
        dict(hypothesis='H4', relationship='Myeloid LGALS3 allocation and β-cell density',
             test='Disease-by-source interaction with donor random slopes', p=float(h4['Pr(>F)'])),
    ])
    summary['p_bh_four_omnibus_tests'] = multipletests(summary.p, method='fdr_bh')[1]
    summary.to_csv(OUT/'primary_disease_comparison_tests.csv',index=False)
    pairs = pd.read_csv(OUT/'h1/donor_pairwise_disease_tests.csv')
    nd_t1d = pairs.loc[pairs.group_1.eq('ND') & pairs.group_2.eq('T1D')].iloc[0]
    figures = {
        'H1': next((OUT/'h1').glob('*.png')),
        'H2': OUT/'h2/coliv_disease_comparison.png',
        'H3': OUT/'h3/h3_density_difference_by_disease.png',
        'H4': OUT/'h4/source_association_by_disease.png',
    }
    legends = {h: next((OUT/h.lower()).glob('figure_legend.*')) for h in figures}
    lines = [
        '# Direct comparisons across disease states — 26 October 2026', '',
        'The primary comparisons are ND versus Aab+, ND versus T1D, and Aab+ versus T1D. '
        'These tests ask whether the specified local relationships differ across disease states. '
        'Pooled associations and tests conducted separately within each state are supporting results; '
        'their significance does not establish a difference between states.', '',
        'The existing corrected cell counts, endocrine objects, marker measurements, ColIV measurements, '
        'and structural geometries were reused. No app data or image pipeline was rebuilt.', '',
        '## Results', '',
        '| Comparison across disease states | Overall p | BH-adjusted p across four overall tests |',
        '|---|---:|---:|',
    ]
    for r in summary.itertuples():
        lines.append(f'| {r.hypothesis}: {r.relationship} | {r.p:.3g} | {r.p_bh_four_omnibus_tests:.3g} |')
    lines += ['',
        f'H1 shows the clearest disease-related pattern: the adjusted immune-allocation association '
        f'is negative in ND and positive in T1D. The direct ND–T1D comparison has '
        f'BH-adjusted p={nd_t1d.p_bh:.3g} across the three H1 pairwise comparisons. '
        f'The H1 overall p={summary.iloc[0].p:.3g} becomes '
        f'{summary.iloc[0].p_bh_four_omnibus_tests:.3g} when adjusted across all four overall tests. '
        'These are different comparison families; the figure asterisk refers to the three pairwise '
        'comparisons, not the four-test adjustment. Neither H1 comparison involving Aab+ is significant.', '',
        'H2, H3, and H4 do not establish differences in their respective relationships across disease '
        'states. This does not establish that the relationships are identical. Earlier positive '
        'within-state or pooled results answer different questions.', '',
        '## Figures and manuscript legends', '',
        'Regression figures use standardized adjusted coefficients (β), not correlation coefficients. '
        'They do not report effects per ten percentage points. A positive β indicates a positive '
        'adjusted relationship and a negative β a negative relationship. The figures omit hypothesis '
        'identifiers and statistical results from titles; significant pairwise comparisons are marked '
        'by brackets and asterisks. Each legend specifies the actual measurements and statistical test.', '',
    ]
    for h in figures:
        legend = legends[h].read_text().strip()
        legend = re.sub(r'^# Figure legend\s*', '', legend)
        lines += [f'### {h}', '', f'![{h} disease comparison]({figures[h]})', '', legend, '']
    lines += [
        '## Supplementary paired donor figures', '',
        'The revised paired figures explicitly compare “Peri only” with “Core present” immune-cell '
        'locations. The presence plot measures the percentage of endocrine objects containing at '
        'least one typed β cell. The density plot measures each donor’s mean of individual object '
        'β-cell densities, with β-negative objects included as zero. These are paired comparisons, '
        'not correlation plots and not the primary disease-state tests.', '',
        f'[Both figures with full manuscript legends]({OLD/"h1/figure_legends.md"}). '
        'PNG, PDF, and SVG versions are saved beside the legends. Their underlying summaries and '
        'paired-test results are unchanged.', '',
        '## Methods and implementation decisions', '',
        'The primary disease comparison must use variation between donors. For H1 and H2, donor '
        'random-intercept/random-slope models reached a singular covariance estimate; the independent '
        'slope versions estimated zero slope variance. Their object-scale tests are retained as '
        'diagnostics. The reported disease comparisons instead use one covariate-adjusted linear '
        'regression coefficient per donor and compare the 20 donor coefficients by Kruskal–Wallis, '
        'followed by all three Wilcoxon rank-sum comparisons with BH adjustment. H3 compares one '
        'observed compartment-density difference per donor using the same group-test family. '
        'Rank-sum tests use exact calculations for untied pairs and tie-corrected approximations '
        'with continuity correction when ties preclude the exact test.', '',
        'H4 retains the nonsingular donor random-intercept/random-slope model. Its disease-by-source '
        'interaction and pairwise slope differences use Satterthwaite degrees of freedom reflecting '
        'donor variation. In H1, H2, and H4, continuous predictors are standardized within donor and '
        'the outcome across all included objects. This separates the local relationship from '
        'differences in average exposure between donors.', '',
        'These method families are represented in the supplied papers: '
        '[Barlow, PDF p22](/data/panc_CODEX/Islet-Explorer/literature/Barlow2025.pdf) '
        '(within-donor predictor standardization, hierarchical linear models and donor random slopes); '
        '[Richardson, PDF p8](/home/smith6jt/Documents/papers/Richardson2023.pdf) '
        '(Kruskal–Wallis disease-group comparisons); '
        '[Bodenmiller, PDF p28](/home/smith6jt/Documents/papers/Bodenmiller2025.pdf) '
        '(donor-aggregated rank-sum comparisons and false-discovery-rate correction). '
        'The specific covariates and exposure definitions remain in each analysis report.', '',
        'These are revised comparisons requested after the earlier analysis, not a preregistered '
        'confirmation. They estimate associations in the existing 20 donors. A useful next step is '
        'to examine the immune lineages contributing to the ND–T1D difference using the same '
        'disease-comparison design.', '',
        '## Files and verification', '',
        'Each hypothesis folder contains the figure, separate legend, model/test tables, and '
        'analysis details. H1/H2 regressions use all 20 donors, H3 uses all 20 paired donor '
        'measurements, and H4 uses 3,681 measurable endocrine objects from 20 donors. Donor 78 '
        'is absent. Source hashes, model diagnostics and reproduced numerical checks are '
        'retained. The original H1 data/model files and previous analysis results remain available.', '',
    ]
    report='\n'.join(lines)
    for path in re.findall(r'\]\((/[^)]+)\)',report):
        assert Path(path).exists(), path
    (OUT/'README.md').write_text(report)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
