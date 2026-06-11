# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
86-dimensional interface featurization (paper Supplementary Methods).
"""

from __future__ import division, print_function

import copy
import json
import sys
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from vaccine_advance_core.featurization.tally_features import (  # noqa: E402
    chemical_class_combinations_reversible,
    size_class_combinations_reversible,
    tally_into_list,
)
from local_pipeline.common.pairs_utils import mutate_pairs_multiple  # noqa: E402
from local_pipeline.common.sequence_utils import diff_seqs  # noqa: E402

DEFAULT_FEATURE_TYPES = [
    chemical_class_combinations_reversible,
    size_class_combinations_reversible,
]


def _feature_values_from_pairs(pairs, mutations):
    mutant_pairs = []
    feature_names = []
    feature_values = []
    for mi in mutations:
        pairs_tmp = copy.deepcopy(pairs)
        pairs_tmp = mutate_pairs_multiple(pairs_tmp, mi)
        feat_vals, feat_names = tally_into_list(pairs_tmp, DEFAULT_FEATURE_TYPES)
        if not feature_names:
            feature_names = feat_names
        feature_values.append(feat_vals)
    return feature_values, feature_names


def load_interface_pairs(path):
    with open(path, 'r') as f:
        data = json.load(f)
    pairs = [tuple(tuple(x) for x in edge) for edge in data['pairs']]
    return pairs, data


def featurize_sequences_df(df, master_seq, interface_pairs, chain_designator=''):
    passthrough_cols = [c for c in df.columns if c != 'sequence']
    records = []
    feature_names = None
    for _, row in df.iterrows():
        muts = diff_seqs(row['sequence'], master_seq, chain_name=chain_designator)
        feat_vals, feat_names = _feature_values_from_pairs(interface_pairs, [muts])
        wt_vals, _ = _feature_values_from_pairs(interface_pairs, [[]])
        combined = feat_vals[0] + wt_vals[0]
        if feature_names is None:
            feature_names = feat_names + ['wt_' + fn for fn in feat_names]
        rec = {
            'sequence': row['sequence'],
            'num_mutations': len(muts),
            'mutation_str': ','.join(
                ['({},{},{},{})'.format(m[0], m[1], m[2], m[3]) for m in muts]
            ),
            'mutationHumanReadable': ','.join(
                ['{}{}{}'.format(m[2], m[1], m[3]) for m in muts]
            ),
        }
        for col in passthrough_cols:
            if col not in rec:
                rec[col] = row[col]
        for fn, fv in zip(feature_names, combined):
            rec[fn] = fv
        records.append(rec)
    return pd.DataFrame(records)
