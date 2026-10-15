#!/usr/bin/env python3
"""Evaluate reference-label transfer with donors held out of training.

This is an explicitly new method, not recovery of the missing scimap rules.
Only marker intensities are predictors. Neither donor identity, diagnosis nor
spatial coordinates are predictors. No new-cell labels or production data are
written. Agreement with historical labels is not biological ground truth.
"""
import argparse
import datetime
import json
from pathlib import Path

import h5py
import numpy as np
import sklearn
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import classification_report, confusion_matrix
from threadpoolctl import threadpool_limits

from prepare_analysis import BASELINE, DEFAULT_ROOT, sha256


def read_vector(group, name):
    value = group[name]
    if isinstance(value, h5py.Group):
        categories = value['categories'].asstr()[:]
        codes = value['codes'][:]
        if (codes < 0).any():
            raise ValueError(f'Missing values in {name}')
        return categories[codes]
    return value.asstr()[:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, default=BASELINE)
    parser.add_argument('--output', type=Path, default=DEFAULT_ROOT / 'reports/phenotype_transfer')
    parser.add_argument('--test-per-donor', type=int, default=5000)
    parser.add_argument('--train-per-donor-label', type=int, default=500)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--feature-space', choices=['raw_log', 'rescaled'], default='raw_log')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260911)
    with h5py.File(args.baseline, 'r') as source:
        donors = read_vector(source['obs'], 'donor_id')
        labels = read_vector(source['obs'], 'phenotype')
        marker_names = source['var/_index'].asstr()[:].tolist()
        # Read the original log1p(raw MFI) layer: the four new raw imports use
        # this identical transform, unlike the baseline's rescaled X matrix.
        values = source['layers/log' if args.feature_space == 'raw_log' else 'X'][:].astype(np.float32)
    if values.shape != (len(labels), len(marker_names)) or not np.isfinite(values).all():
        raise ValueError('Invalid reference marker matrix')
    classes = np.unique(labels)
    donor_ids = sorted(np.unique(donors))
    groups = [donor_ids[i::5] for i in range(5)]
    sample_train, sample_test = {}, {}
    for donor in donor_ids:
        members = np.flatnonzero(donors == donor)
        sample_test[donor] = rng.choice(members, size=min(args.test_per_donor, len(members)), replace=False)
        sample_train[donor] = np.concatenate([
            rng.choice(indices, size=min(args.train_per_donor_label, len(indices)), replace=False)
            for label in classes if len(indices := members[labels[members] == label])])
    params = dict(max_iter=100, max_leaf_nodes=31, max_depth=6, learning_rate=0.1,
                  l2_regularization=1., early_stopping=False, random_state=20260911)
    # Fix screening criteria before scoring. These are technical criteria for
    # deciding whether to investigate transfer further, not scientific approval.
    criteria = {'accuracy_at_least': 0.97, 'macro_f1_at_least': 0.90,
                'recall_at_least_for_labels_with_50_test_cells': 0.85}
    settings = {'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'method': 'HistGradientBoostingClassifier reference-label transfer, five donor-disjoint folds',
                'parameters': params, 'feature_space': args.feature_space,
                'marker_names': marker_names, 'test_per_donor': args.test_per_donor,
                'train_per_donor_label': args.train_per_donor_label, 'screening_criteria': criteria,
                'fold_donors': groups, 'script_sha256': sha256(__file__), 'sklearn_version': sklearn.__version__,
                'baseline_path': str(args.baseline), 'baseline_bytes': args.baseline.stat().st_size}
    (args.output / 'design.json').write_text(json.dumps(settings, indent=2) + '\n')
    observed, predicted, scores, folds = [], [], [], []
    for i, test_donors in enumerate(groups):
        train_donors = sorted(set(donor_ids) - set(test_donors))
        train = np.concatenate([sample_train[donor] for donor in train_donors])
        test = np.concatenate([sample_test[donor] for donor in test_donors])
        if set(donors[train]) & set(donors[test]):
            raise ValueError('Donor leakage between training and evaluation')
        print(f'Fold {i+1}: {len(train):,} training cells; held-out donors {test_donors}', flush=True)
        with threadpool_limits(limits=args.threads):
            model = HistGradientBoostingClassifier(**params).fit(values[train], labels[train])
            probabilities = model.predict_proba(values[test])
        prediction = model.classes_[probabilities.argmax(axis=1)]
        observed.extend(labels[test].tolist())
        predicted.extend(prediction.tolist())
        scores.extend(probabilities.max(axis=1).tolist())
        fold = {'fold': i+1, 'train_donors': train_donors, 'test_donors': test_donors,
                'training_cells': len(train), 'test_cells': len(test),
                'classification': classification_report(labels[test], prediction, labels=classes,
                                                          output_dict=True, zero_division=0)}
        folds.append(fold)
        (args.output / f'fold_{i+1}.json').write_text(json.dumps(fold, indent=2) + '\n')
        print(f"Fold {i+1}: accuracy {np.mean(prediction == labels[test]):.4f}", flush=True)
    report = classification_report(observed, predicted, labels=classes, output_dict=True, zero_division=0)
    accuracy = float(np.mean(np.asarray(observed) == np.asarray(predicted)))
    passing = (accuracy >= criteria['accuracy_at_least'] and report['macro avg']['f1-score'] >= criteria['macro_f1_at_least']
               and all(report[label]['recall'] >= criteria['recall_at_least_for_labels_with_50_test_cells']
                       for label in classes if report[label]['support'] >= 50))
    high = np.asarray(scores) >= 0.9
    summary = {'design': settings, 'classification': report, 'accuracy': accuracy,
               'confusion_labels': classes.tolist(), 'confusion_matrix': confusion_matrix(observed, predicted, labels=classes).tolist(),
               'screening_passed': bool(passing), 'status': 'eligible_for_further_review' if passing else 'rejected_for_automatic_integration',
               'score_ge_0_9_fraction': float(high.mean()),
               'score_ge_0_9_accuracy': float(np.mean(np.asarray(observed)[high] == np.asarray(predicted)[high])) if high.any() else None,
               'new_cell_labels_written': False, 'original_method_recovered': False,
               'limits': ['Reference agreement does not validate the biology of historical labels.',
                          'Classifier scores are not calibrated probabilities.',
                          'Saved rescaled features require independent reproduction of the original transform before new-cell use.',
                          'Original and new acquisitions also differ in segmentation version and require image review.',
                          'No thresholds were relaxed using these evaluation results.']}
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'accuracy': accuracy, 'macro_f1': report['macro avg']['f1-score'], 'screening_passed': bool(passing)}), flush=True)


if __name__ == '__main__':
    main()
