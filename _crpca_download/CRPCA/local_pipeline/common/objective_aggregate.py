# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Aggregate multipoint objectives (Step 07)."""

from __future__ import division, print_function

import ast

import numpy as np
import pandas as pd

from local_pipeline.common.scoring import SCORE_COLUMNS


def build_single_point_lookup(scores_df):
    lookup = {}
    for _, r in scores_df.iterrows():
        key = '{}{}{}'.format(r['original_aa'], int(r['location']), r['mutant_aa'])
        lookup[key] = {c: r[c] for c in SCORE_COLUMNS}
    return lookup


def parse_mutations_from_row(row):
    if 'mutationHumanReadable' in row and isinstance(row['mutationHumanReadable'], str):
        return [p.strip() for p in row['mutationHumanReadable'].split(',') if p.strip()]
    if 'mutations' in row and isinstance(row['mutations'], str):
        muts = ast.literal_eval(row['mutations'])
        return ['{}{}{}'.format(m[2], m[1], m[3]) for m in muts]
    raise ValueError('无法解析突变')


def aggregate_objectives(sequences_df, single_point_scores_df):
    lookup = build_single_point_lookup(single_point_scores_df)
    rows = []
    for _, row in sequences_df.iterrows():
        keys = parse_mutations_from_row(row)
        agg = {c: 0.0 for c in SCORE_COLUMNS}
        # --- PATCH: 按工具独立判定有效性。某工具在某个构成突变上缺失（NaN）时，
        # 仅令该工具的 sum 为 NaN，不再连累其他工具（支持 sfe 等留空的情况）。
        tool_valid = {c: True for c in SCORE_COLUMNS}
        found = True
        for k in keys:
            if k not in lookup:
                found = False
                break
            for c in SCORE_COLUMNS:
                v = lookup[k][c]
                if pd.isna(v):
                    tool_valid[c] = False
                else:
                    agg[c] += float(v)
        rec = {
            'sequence': row['sequence'],
            'num_mutations': len(keys),
            'mutationHumanReadable': ','.join(keys),
            'all_single_points_found': found,
        }
        if 'pred_mean' in row:
            rec['pred_mean'] = row['pred_mean']
        if 'mei_score' in row:
            rec['mei_score'] = row['mei_score']
        if 'bo_round' in row:
            rec['bo_round'] = row['bo_round']
        for c in SCORE_COLUMNS:
            rec['sum_' + c] = agg[c] if (found and tool_valid[c]) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)
