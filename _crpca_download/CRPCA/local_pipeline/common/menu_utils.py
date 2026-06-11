# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""从单点打分 CSV 与 master 序列构建 allowed_mutations 菜单。"""

from __future__ import division, print_function

import json


def build_allowed_mutations_from_scores(scores_df, master_seq):
    """
    :return: list of [location, original_aa, [mutant_aa, ...]]
    """
    df = scores_df.copy()
    df['location'] = df['location'].astype(int)
    allowed = []
    for loc, grp in df.groupby('location'):
        orig_vals = grp['original_aa'].unique()
        if len(orig_vals) != 1:
            raise ValueError('位点 {} 存在多个 original_aa: {}'.format(loc, orig_vals))
        orig = orig_vals[0]
        if loc < 1 or loc > len(master_seq):
            raise ValueError('位点 {} 超出 master 长度 {}'.format(loc, len(master_seq)))
        master_aa = master_seq[loc - 1]
        if master_aa != orig:
            raise ValueError(
                '位点 {}: CSV original_aa={} 与 master 上 {} 不一致'.format(
                    loc, orig, master_aa)
            )
        mutants = sorted(set(grp['mutant_aa'].astype(str).tolist()))
        allowed.append([int(loc), str(orig), mutants])
    allowed.sort(key=lambda x: x[0])
    return allowed


def save_allowed_mutations(allowed, path):
    with open(path, 'w') as f:
        json.dump(allowed, f, indent=2)
