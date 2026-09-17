"""Apply validated retained rules and separately audit archival assignments.

Satisfying a threshold is a computational consistency check, not independent
biological validation. Unvalidated dependencies cannot certify a new label.
"""
import csv
import numpy as np
from common import RETAINED

ALIASES = {'SMA+':'Smooth Muscle','Blood Vessel':'Blood vessel',
           'CD8a Tcell':'CD8 T cell','CD4 Tcell':'CD4 T cell'}

def read_rules(path):
    with open(path) as stream: rows = list(csv.reader(stream))
    markers = rows[0][2:]; rules = {}
    for row in rows[1:]:
        if len(row)<2 or not row[1]: continue
        rules[row[1]] = (row[0], {m:v for m,v in zip(markers,row[2:]) if v})
    return rules

def dependencies(name, rules, seen=()):
    name = ALIASES.get(name,name)
    if name in seen: raise ValueError('Cyclic phenotype parent rules')
    parent, tests = rules[name]
    chain = [] if parent == 'all' else dependencies(parent,rules,seen+(name,))
    return chain + [tests]

def rule_mask(log_values, markers, chain, gates):
    keep = np.isfinite(log_values).all(axis=1)
    for tests in chain:
        anypos = []
        for marker,test in tests.items():
            if marker not in gates or not np.isfinite(gates[marker]):
                return np.zeros(len(log_values),bool)
            positive = log_values[:,markers.index(marker)] > gates[marker]
            if test == 'anypos': anypos.append(positive)
            elif test == 'pos': keep &= positive
            elif test == 'allneg': keep &= ~positive
            else: raise ValueError(f'Unknown gate operator {test}')
        if anypos: keep &= np.logical_or.reduce(anypos)
    return keep

def audit(labels, logs, markers, rules, gates, validated_markers):
    records = []; masks = {}
    for label in RETAINED:
        chain = dependencies(label,rules)
        deps = sorted(set().union(*(set(x) for x in chain)))
        matches = rule_mask(logs,markers,chain,gates)
        existing = labels == label
        computable = all(m in gates and np.isfinite(gates[m]) for m in deps)
        missing_thresholds = {m for m in deps if m not in gates or not np.isfinite(gates[m])}
        missing = sorted((set(deps)-set(validated_markers)) | missing_thresholds)
        # Ductal/Acinar share a negative-only rule in the supplied CSV.
        specific = label != 'Ductal' or any(
            tests.get(m) in {'pos', 'anypos'} for tests in chain for m in ['CK19', 'panCK'])
        validated = not missing and specific
        # New donors do not have archival labels. Rule application must not
        # depend on a previous assignment; keep the archival comparison below.
        masks[label] = matches if validated else np.zeros(len(labels),bool)
        records.append(dict(phenotype=label, existing_count=int(existing.sum()),
                            gate_consistent_count=int((existing & matches).sum()) if computable else np.nan,
                            rule_positive_count=int(matches.sum()) if computable else np.nan, validated=validated,
                            rule_computable=computable,
                            assigned_count=int(masks[label].sum()),
                            missing_validated_gates=';'.join(missing),
                            reason='' if validated else ('ductal_rule_lacks_CK19_panCK_specificity' if not specific
                                                       else 'upstream_gate_validation_unavailable')))
    return records, masks
