#!/usr/bin/env python3
"""Numerical H6 reporting and manuscript legends; no branch assignments."""
from pathlib import Path
import pandas as pd

OUT=Path('output/h5_h10_original_plan_20261026/h6')

def main():
    basis=('The dataset comprises 721 eligible endocrine objects (371 ND, 314 Aab+, 36 T1D) from 19 donors '
           '(7 ND, 7 Aab+, 5 T1D). Donor 6563 has no eligible graph and remains in the coverage table with unavailable values. '
           'Eligibility requires at least ten typed endocrine cells, at least two β cells, at least two α and δ cells combined, '
           'and at least one graph edge. The original current cell labels and core memberships are retained. '
           'An undirected edge connects endocrine centroids within 20 µm. For each object, 1,999 random label permutations '
           'preserve its graph geometry and every cell-type count (master seed 20261026 with stable object-specific streams). '
           'Combined adjacency is the number of β–α plus β–δ edges divided by all endocrine–endocrine edges. '
           'Excess adjacency is the observed fraction minus its permutation mean. All 721 graphs have nonzero combined null variance. ')
    adjustment=('The fixed covariates are β-cell fraction, δ-cell fraction and other-endocrine fraction with α cells as reference, '
                'log endocrine-core area, log1p nearby immune-cell count and log1p core immune-cell count. '
                'Nearby immune cells comprise the core plus its 20-µm peri-object region. Continuous predictors are standardized '
                'within donor and the outcome across included objects. The outcome is mean existing scaled cytoplasmic INS '
                'among retained typed β cells. Coefficients are standardized regression coefficients, not correlations. ')
    (OUT/'observed_vs_permuted_adjacency_legend.txt').write_text(
        'Observed and permuted endocrine adjacency. '+basis+
        'Each connected pair is one donor’s arithmetic mean across eligible objects; the left point is observed adjacency '
        'and the right point is the permutation-mean adjacency, both as percentages of all endocrine edges. '
        'Exact two-sided paired Wilcoxon signed-rank p values are 0.015625 for ND, 0.015625 for Aab+, and 0.812500 for T1D. '
        'BH adjustment across these three disease-specific paired comparisons gives 0.023438, 0.023438 and 0.812500. '
        'These tests compare adjacency with relabeling within the fixed graphs; the conditional INS regressions are reported separately.\n')
    (OUT/'excess_adjacency_by_disease_legend.txt').write_text(
        'Supplementary donor comparison of excess endocrine adjacency. '+basis+
        'Each point is one donor’s mean excess fraction, expressed in percentage points. Boxes show the median and interquartile range; '
        'whiskers reach observations within 1.5 interquartile ranges, and every donor is displayed. The horizontal dashed line marks zero. '
        'Median donor values are −3.8317 percentage points for ND, −4.8722 for Aab+, and 0.4729 for T1D. '
        'The Kruskal–Wallis statistic is 3.475489, with p=0.175917. Two-sided exact Wilcoxon rank-sum p values for ND–Aab+, ND–T1D '
        'and Aab+–T1D are 0.455711, 0.148990 and 0.148990; BH adjustment across these three pairs gives 0.455711, 0.223485 and 0.223485. '
        'This supplementary measurement comparison was calculated after the original conditional-model specification and is not the original INS-association test.\n')
    (OUT/'adjacency_insulin_association_legend.txt').write_text(
        'Adjusted endocrine adjacency and retained β-cell insulin. '+basis+adjustment+
        '\nLeft: each point is a full-rank donor-specific regression coefficient with positive residual degrees of freedom. '
        'There are seven ND, seven Aab+, and one T1D coefficients. Boxes show median and interquartile range, whiskers 1.5 interquartile ranges, '
        'and all points are shown; no box is drawn for the single T1D coefficient. Donor 6533 contributes the T1D coefficient. '
        'Donors 6520, 6534, 6550 and 6551 have respectively 7, 6, 3 and 4 eligible objects for eight requested model columns, '
        'and zero residual degrees of freedom. No nonconstant covariates were removed. The all-donor mean coefficient is 0.087428 '
        '(Student-t 95% CI −0.044399 to 0.219255), and the median is 0.011299; the two-sided signed-rank test against zero gives p=0.151428. '
        'ND–Aab+ rank-sum p=0.259324 (BH p=0.777972 using the planned three-pair family). The T1D pairwise tests and three-state global test '
        'are not reported because there is only one independently estimable T1D coefficient. '
        '\nRight: points and bars are disease-specific slopes and 95% CIs from the supporting random-intercept model with the same fixed covariates. '
        'The initial correlated random-slope model has intercept–slope correlation −1; the independent random-slope model has contact-slope variance zero. '
        'The displayed random-intercept model therefore assumes no residual donor slope heterogeneity within disease states. It uses all 721 objects '
        'and 19 donors. Slopes are −0.009417 (95% CI −0.065614 to 0.046780) for ND, 0.019181 (−0.041348 to 0.079710) for Aab+, and '
        '0.019042 (−0.170126 to 0.208210) for T1D. The contact-by-disease F test gives p=0.776851; BH adjustment of the three slope differences gives '
        '0.998893 for every pair. In the separate common-slope model, β=0.004416 (95% CI −0.036745 to 0.045578), p=0.833221. '
        'Full singular models, supporting models, donor fits and diagnostics are retained.\n')
    sensitivity=pd.read_csv(OUT/'plot_data_sensitivity_estimates.csv')
    lines=[]
    for _,r in sensitivity.iterrows():
        lines.append(f'{r.label}: β={r.standardized_beta:.6f}, 95% CI {r.lower:.6f} to {r.upper:.6f}, p={r.p_value:.6f}; {r.objects} objects and {r.donors} donors.')
    (OUT/'adjacency_analysis_variants_legend.txt').write_text(
        'Conditional coefficients across the original analysis and same-data variants. '+basis+adjustment+
        'Points are common-slope estimates and bars are 95% CIs. A correlated donor intercept/slope model is used where nonsingular; '
        'the independent covariance version is retained if the correlated fit is singular. A supporting random-intercept-only estimate is displayed '
        'where the independent contact-slope variance reaches zero. This selection rule uses covariance diagnostics, not p values. '
        'All models and their covariance variants are supplied. Separate α and δ models require at least two cells of the respective partner type '
        'and nonzero corresponding permutation variance. The β-density outcome uses typed β cells per endocrine-core mm². '
        'The contact-by-immune row is the interaction between standardized excess adjacency and standardized log1p nearby immune count. '
        'Raw and RedSEA-corrected rows use existing cytoplasmic INS measurements; stronger hormone margin retains typed β cells with INS at least '
        'four times the maximum GCG/SST signal. Nerve/vessel adjustment uses INS-independent B3TUBB and CD31 area fractions in tissue-supported '
        '20-µm peri-object regions. The proximity variant adds the fraction of B3TUBB nerve area within 20 µm of CD31, leaving absent-nerve proximity '
        'unavailable. True-boundary contacts use intersecting existing cell polygons, without dilation or geometry repair; 39,438 edges occur across '
        'the 721 graphs and receive the same count-preserving permutation procedure. Each variant standardizes the outcome over its included objects. '
        'The reported p values are unadjusted individual sensitivity-model values; no significance annotations or cross-variant decisions are displayed.\n\n'+
        '\n'.join(lines)+'\n')
    (OUT/'matched_composition_size_insulin_legend.txt').write_text(
        'Within-stratum INS measurements. Eligible endocrine objects are grouped using pooled tertiles of β-cell fraction, δ-cell fraction and '
        'log endocrine-core area, calculated without INS outcomes. Within each donor and three-variable stratum, both positive and nonpositive '
        'combined excess-adjacency objects are required. Each connected pair shows that donor’s mean of stratum-level mean β-cell INS, '
        'with equal weight for each qualifying stratum. There are 76 qualifying strata and 403 objects: five ND donors (33 strata, 129 nonpositive '
        'and 75 positive objects), six Aab+ donors (39 strata, 106 nonpositive and 80 positive objects), and four T1D donors (4 strata, '
        '7 nonpositive and 6 positive objects). Tertile edges, each object’s assignment, stratum counts and donor values are exported. '
        'These are coarse composition/size strata, not exact composition matches. No additional immune matching, significance tests or multiplicity '
        'adjustment are applied to this descriptive display; immune-adjusted continuous regressions are reported separately. INS units are the '
        'existing scaled cytoplasmic cell signal.\n')
    # Keep the earlier observed-only display, but remove the now-outdated claim
    # that no permutations have been run anywhere in the output directory.
    legacy=OUT/'observed_beta_alpha_delta_adjacency_legend.txt'
    if legacy.exists():
        text=legacy.read_text()
        text=text.replace('No label permutations, excess-over-random inference, disease tests, or multiple-comparison adjustments were performed. Differences can reflect endocrine composition and do not establish the original spatial-null hypothesis.',
            'This observed-only figure contains no statistical tests. The subsequently approved label-permutation measurements and conditional regressions are provided in separate figures and tables.')
        legacy.write_text(text)
    table='\n'.join(f'| {r.label} | {r.standardized_beta:.6f} | {r.lower:.6f}, {r.upper:.6f} | {r.p_value:.6f} | {r.objects} |' for _,r in sensitivity.iterrows())
    readme=f'''# H6 numerical outputs and original-plan execution

This folder contains the approved fixed-geometry, count-preserving label permutations, conditional regressions, original same-data variants, figures and machine-readable outputs. No biological branch is assigned here.

## Dataset and measurement

{basis}

The source polygon sensitivity contains 39,438 contact edges; all 721 objects have edges and nonzero combined permutation variance. No source polygons are invalid. Complete permutation draws are stored in `permutation_null_draws.npz` and `boundary_permutation_null_draws.npz`; each array retains object identifiers and the three statistics. The primary excess uses the permutation mean; the exact combinatorial expected fraction and Monte Carlo mean error are retained for verification.

## Conditional coefficients and disease comparisons

{adjustment}

The original correlated random-slope fit has intercept–slope correlation −1. Its independent covariance version has contact-slope variance zero. The separate donor regressions are full rank with positive residual degrees of freedom in 15 donors: seven ND, seven Aab+, one T1D (6533). Their mean coefficient is 0.087428 (95% Student-t CI −0.044399 to 0.219255); the median is 0.011299 and the signed-rank p value against zero is 0.151428. There are ten positive and five negative coefficients. The ND–Aab+ rank-sum p value is 0.259324 (BH p=0.777972 with the planned three-pair family). No independently estimated T1D pairwise or three-state global test is reported because the T1D coefficient count is one. Donors 6520/6534/6550/6551 have 7/6/3/4 observations for eight requested model columns and zero residual degrees of freedom. Donor 6563 has no eligible graph. Required nonconstant covariates are retained.

The supporting random-intercept model uses all 721 objects and 19 donors, assuming no residual donor slope heterogeneity within disease states. Its common coefficient is 0.004416 (95% CI −0.036745 to 0.045578), p=0.833221. The contact-by-disease interaction p=0.776851; the three BH-adjusted slope-contrast p values are all 0.998893. Its disease-specific slopes and intervals are exported in `plot_data_primary_mixed_model_slopes.csv`. The random-intercept fit does not replace the independent donor-coefficient results. Both the fixed original specification and the covariance-boundary amendment are saved.

## Numerical sensitivity results

Individual sensitivity p values below are unadjusted; no cross-variant branch is assigned. All displayed model fits are nonsingular. Each variant’s earlier singular covariance forms remain saved.

| Analysis | Standardized β | 95% CI | p | Objects |
|---|---:|---|---:|---:|
{table}

## Additional measurements

The observed-versus-permuted donor paired tests give raw p=0.015625/0.015625/0.812500 and BH p=0.023438/0.023438/0.812500 for ND/Aab+/T1D. The three-test family is the disease-specific paired comparisons. The additional across-state excess-adjacency comparison was calculated after the original conditional specification and is retained as supplementary: Kruskal–Wallis p=0.175917; ND–Aab+, ND–T1D and Aab+–T1D rank-sum BH p=0.455711/0.223485/0.223485. It is not the original INS-association test.

The composition/size display includes 76 within-donor strata containing both positive and nonpositive excess-adjacency objects, comprising 403 objects and 15 donors. It uses pooled tertiles of β fraction, δ fraction and log core area; exact bin boundaries and all assignments are saved. This display has no additional hypothesis test.

`original_plan_execution.csv` maps each original primary/secondary/follow-up operation to its output. The original “if supported” and “if a meaningful effect is excluded” conditions are listed without selecting either branch. No meaningful-effect threshold was supplied and no equivalence analysis was executed.

## Reproduction and files

The existing H5 object features and saved endocrine graph inputs are reused. The H7 structural covariates use INS-independent B3TUBB/CD31 measurements. The original cell polygons come from `data/development/uniform-20donor-20261015/cell_geojson/`; no segmentation, typing or pipeline changes were made. The 20-µm graph construction has local precedent in `~/Documents/papers/human_panc_atlas.pdf`, PDF p30. Standardization and donor mixed-model methods follow repository `literature/Barlow2025.pdf`, PDF p22, and `~/Documents/papers/Bodenmiller2025.pdf`, PDF p28. The original label-permutation operation was explicitly approved by the user for this continuation.

```sh
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/research/h6_permutation_measurements_20261026.py
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/research/h6_boundary_contacts_20261026.py
Rscript scripts/research/h6_models_20261026.R
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/research/h6_figures_20261026.py
/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/research/h6_reports_20261026.py
```

Every new figure has PNG/PDF/SVG outputs and a separate `<stem>_legend.txt`. Full object features, permutation measurements, donor-level values, model coefficients, covariance diagnostics, model RDS files and plotted data are retained. `README_descriptive_stage.md`, `descriptive_stage_validation.json` and `descriptive_stage_manifest.json` preserve the earlier descriptive-only stage; the current validation and manifest apply to the completed computations.
'''
    (OUT/'README.md').write_text(readme)
    operations=[
        ('Primary','Count/geometry-preserving label permutations','Executed','permutation_measurements.csv'),
        ('Primary','Combined contact statistic','Executed','h6_object_features.csv'),
        ('Primary','Adjusted main INS association','Executed','donor_direction_tests.csv; model_coefficients.csv'),
        ('Presentation','Disease comparisons of adjusted association','Mixed models executed; independent T1D donor count is one','disease_slope_comparisons.csv; donor_disease_pairwise_tests.csv'),
        ('Primary','Composition/size strata','Coarse within-donor tertile strata displayed','matched_composition_size_strata.csv'),
        ('Secondary','Beta density outcome','Executed','beta_density_analysis_objects.csv; model_coefficients.csv'),
        ('Secondary','Contact-by-immune interaction','Executed','contact_by_immune_analysis_objects.csv; model_coefficients.csv'),
        ('Same-data follow-up','True cell-boundary contacts','Executed from existing cell polygons','boundary_contact_measurements.csv'),
        ('Same-data follow-up','Separate Alpha and Delta','Executed','alpha_analysis_objects.csv; delta_analysis_objects.csv'),
        ('Same-data follow-up','Nerve/vessel conditioning','Density and available H7 proximity covariates executed','nerve_vessel_proximity_adjustment_analysis_objects.csv'),
        ('Condition in original plan','If supported','No branch assigned',''),
        ('Condition in original plan','If meaningful effect excluded','No branch assigned; no equivalence bound supplied or analysis executed',''),
    ]
    pd.DataFrame(operations,columns=['plan_component','operation_or_condition','execution_record','output']).to_csv(OUT/'original_plan_execution.csv',index=False)

if __name__=='__main__':main()
