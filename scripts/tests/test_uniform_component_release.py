"""Metadata contracts only; no synthetic typing or image evidence."""
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import assemble_uniform as assembly
from composition_embedding import endocrine_labels


def test_resolved_group_uses_the_actual_four_name_registry_without_status_inference():
    lab = assembly.ID.labels_table(assembly.RG.load())
    names = list(assembly.NAMED_ENDOCRINE)
    group = assembly.endocrine_group(lab, names + ['Unclassified', 'Exocrine'], {},
                                    assembly.RESOLVED_COMPONENT_SCHEMA)
    assert set(group) == set(names)
    assert len(group) == 4
    assert 'Unclassified' not in group
    for missing in names:
        with pytest.raises(ValueError, match='all four'):
            assembly.endocrine_group(lab, [x for x in names if x != missing], {},
                                     assembly.RESOLVED_COMPONENT_SCHEMA)
    for unresolved in ('Islet', 'EndocrineCluster', 'EndocrineSingle', 'Endocrine'):
        with pytest.raises(ValueError, match='unnamed'):
            assembly.endocrine_group(lab, names + [unresolved], {}, assembly.RESOLVED_COMPONENT_SCHEMA)


def test_historical_group_is_preserved():
    lab = assembly.ID.labels_table(assembly.RG.load())
    names = list(assembly.NAMED_ENDOCRINE)
    comp = {label: {('Epithelial', 'Endocrine')} for label in names}
    labels = names + list(assembly.ENDOCRINE_STARTS) + ['Islet', 'Tissue']
    assert assembly.endocrine_group(lab, labels, comp, 'uniform.component_selection/1') == labels[:-2]


def test_resolved_embedding_rejects_an_unnamed_or_missing_denominator_member():
    names = list(assembly.NAMED_ENDOCRINE)
    source = SimpleNamespace(uns={'endocrine_group_contract': 'named_four/1', 'endocrine_labels': names})
    assert endocrine_labels(source) == names
    for labels in (names + ['Unclassified'], names[:-1], names[:-1] + ['EndocrineSingle']):
        source.uns['endocrine_labels'] = labels
        with pytest.raises(ValueError, match='exactly the four'):
            endocrine_labels(source)
