# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Structural Fluctuation Estimation (SFE) — 从 Rosetta Flex ddG ensemble 聚合最终分数。

论文 Supplementary Methods:
  1. 对 forward（WT 构象）与 reverse（突变构象）两组 Flex ddG 分别去 outlier
  2. 对四分位范围内结果取平均
  3. SFE = (ddG_forward - ddG_reverse) / 2

输入：宽表 CSV，mutation + chain + forward_01..forward_21 + reverse_01..reverse_21（42 个 ddG）。
"""

from __future__ import division, print_function

import numpy as np
import pandas as pd

DEFAULT_N_CONFORMATIONS = 21


def interquartile_mean(values, iqr_multiplier=1.5, mode='tukey'):
    """
    对一组 ddG 做 outlier 过滤后取平均。

    :param mode:
        'tukey' — 去掉 Q1-1.5*IQR 与 Q3+1.5*IQR 之外的值，对剩余取均值（默认）
        'strict-iqr' — 仅保留 [Q1, Q3] 内的值再取均值
    """
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        raise ValueError('ddG 数组不能为空')
    if np.any(~np.isfinite(arr)):
        raise ValueError('ddG 数组含 NaN/Inf')

    q1, q3 = np.percentile(arr, [25.0, 75.0])
    iqr = q3 - q1
    lower = q1 - iqr_multiplier * iqr
    upper = q3 + iqr_multiplier * iqr

    if mode == 'tukey':
        mask = (arr >= lower) & (arr <= upper)
    elif mode == 'strict-iqr':
        mask = (arr >= q1) & (arr <= q3)
    else:
        raise ValueError('未知 mode: {}'.format(mode))

    used = arr[mask]
    if used.size == 0:
        used = arr

    return {
        'mean': float(np.mean(used)),
        'n_total': int(arr.size),
        'n_used': int(used.size),
        'q1': float(q1),
        'q3': float(q3),
        'iqr_lower': float(lower),
        'iqr_upper': float(upper),
    }


def compute_sfe_ddg(
        forward_ddgs,
        reverse_ddgs,
        iqr_multiplier=1.5,
        mode='tukey',
        n_conformations=None,
):
    """计算单个突变的 SFE ddG。"""
    fwd = np.asarray(forward_ddgs, dtype=float)
    rev = np.asarray(reverse_ddgs, dtype=float)

    if n_conformations is not None:
        if fwd.size != n_conformations or rev.size != n_conformations:
            raise ValueError(
                '期望每组 {} 个 ddG，实际 forward={}, reverse={}'.format(
                    n_conformations, fwd.size, rev.size))

    fwd_stats = interquartile_mean(fwd, iqr_multiplier=iqr_multiplier, mode=mode)
    rev_stats = interquartile_mean(rev, iqr_multiplier=iqr_multiplier, mode=mode)

    sfe_ddg = (fwd_stats['mean'] - rev_stats['mean']) / 2.0

    return {
        'sfe_ddg': float(sfe_ddg),
        'forward_mean': fwd_stats['mean'],
        'reverse_mean': rev_stats['mean'],
        'forward_n_used': fwd_stats['n_used'],
        'reverse_n_used': rev_stats['n_used'],
        'forward_n_total': fwd_stats['n_total'],
        'reverse_n_total': rev_stats['n_total'],
        'forward_q1': fwd_stats['q1'],
        'forward_q3': fwd_stats['q3'],
        'reverse_q1': rev_stats['q1'],
        'reverse_q3': rev_stats['q3'],
    }


def _wide_ddg_column_names(n_conformations):
    fwd = ['forward_{:02d}'.format(i) for i in range(1, n_conformations + 1)]
    rev = ['reverse_{:02d}'.format(i) for i in range(1, n_conformations + 1)]
    return fwd, rev


def _empty_string_mask(series):
    return series.isna() | (series.astype(str).str.strip() == '')


def load_flex_ddg_wide(csv_path, n_conformations=DEFAULT_N_CONFORMATIONS):
    """
    宽表 CSV：mutation, chain, forward_01..forward_NN, reverse_01..reverse_NN
    """
    df = pd.read_csv(csv_path)
    for col in ('mutation', 'chain'):
        if col not in df.columns:
            raise ValueError('宽表须含 {} 列'.format(col))

    empty_chain = _empty_string_mask(df['chain'])
    if empty_chain.any():
        raise ValueError(
            'chain 列不能为空（{} 行缺失）'.format(int(empty_chain.sum())))

    fwd_cols, rev_cols = _wide_ddg_column_names(n_conformations)
    missing = [c for c in fwd_cols + rev_cols if c not in df.columns]
    if missing:
        raise ValueError('宽表缺少列: {}'.format(missing))

    rows = []
    for _, row in df.iterrows():
        fwd = row[fwd_cols].astype(float).tolist()
        rev = row[rev_cols].astype(float).tolist()
        rows.append({
            'mutation': row['mutation'],
            'chain': str(row['chain']).strip(),
            'forward_ddgs': fwd,
            'reverse_ddgs': rev,
        })
    return rows


def aggregate_sfe_table(
        records,
        iqr_multiplier=1.5,
        mode='tukey',
        n_conformations=DEFAULT_N_CONFORMATIONS,
        invert_for_sampling=False,
):
    """批量计算 SFE，返回 DataFrame。"""
    out_rows = []
    for rec in records:
        stats = compute_sfe_ddg(
            rec['forward_ddgs'],
            rec['reverse_ddgs'],
            iqr_multiplier=iqr_multiplier,
            mode=mode,
            n_conformations=n_conformations,
        )
        row = {'mutation': rec['mutation'], 'chain': rec['chain'], **stats}
        row['sfe'] = -row['sfe_ddg'] if invert_for_sampling else row['sfe_ddg']
        out_rows.append(row)
    return pd.DataFrame(out_rows)
