# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Score transforms for the local antibody design pipeline."""

from __future__ import division, print_function

import numpy as np
import pandas as pd

# Paper Supplementary Methods: generalized logistic for mutation sampling.
LOGISTIC_A = 1000.0
LOGISTIC_B = 5.0
LOGISTIC_C = 2.0

SCORE_COLUMNS = ['sfe', 'fep', 'rosetta_flex', 'foldx', 'abnativ']
TOOL_COLUMNS = SCORE_COLUMNS  # alias

SINGLE_POINT_REQUIRED_COLUMNS = (
    ['mutation', 'chain', 'location', 'original_aa', 'mutant_aa'] + SCORE_COLUMNS
)


def _empty_string_mask(series):
    return series.isna() | (series.astype(str).str.strip() == '')


def validate_single_point_scores_df(df, antibody_chains=None):
    """Validate single-point score CSV columns and chain values."""
    missing = [c for c in SINGLE_POINT_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError('单点打分 CSV 缺少列: {}'.format(missing))

    empty_mutation = _empty_string_mask(df['mutation'])
    if empty_mutation.any():
        raise ValueError(
            'mutation 列不能为空（{} 行缺失）'.format(int(empty_mutation.sum())))

    empty_chain = _empty_string_mask(df['chain'])
    if empty_chain.any():
        raise ValueError(
            'chain 列不能为空（{} 行缺失）'.format(int(empty_chain.sum())))

    if antibody_chains is not None:
        ab_set = set(antibody_chains)
        bad = ~df['chain'].astype(str).str.strip().isin(ab_set)
        if bad.any():
            bad_chains = sorted(df.loc[bad, 'chain'].astype(str).str.strip().unique())
            raise ValueError(
                'chain 须属于抗体链 {}；发现: {}'.format(
                    list(antibody_chains), bad_chains))


def logistic_transform(score):
    """Map a raw tool score to an unnormalized probability mass."""
    score = np.asarray(score, dtype=float)
    denom = np.power(1.0 + LOGISTIC_A * np.exp(-score * LOGISTIC_B), 1.0 / LOGISTIC_C)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = 1.0 / denom
    return out


def combine_tool_scores(row, tool_columns=None):
    """Sum logistic-transformed scores across tools for one mutation."""
    if tool_columns is None:
        tool_columns = TOOL_COLUMNS
    total = 0.0
    for col in tool_columns:
        val = row[col]
        if pd.isna(val):
            continue
        total += float(logistic_transform(val))
    return total


def normalize_sampling_weights(df, weight_column='samplingWeight'):
    """Normalize weights to sum to 1."""
    total = df[weight_column].sum()
    if total <= 0:
        raise ValueError('Total sampling weight is non-positive; check score orientation.')
    df = df.copy()
    df[weight_column] = df[weight_column] / total
    return df
