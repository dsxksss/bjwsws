#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 02 — 从单点突变五工具打分 CSV 计算序列生成器采样权重。

对应论文 2.2（已完成单点扫描）→ 为 2.3 Step B 准备 samplingWeight。

用法:
  python step02_compute_sampling_weights.py \\
    --input scores/single_point_scores.csv \\
    --output scores/sampling_weights.csv
"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.scoring import (  # noqa: E402
    SCORE_COLUMNS,
    combine_tool_scores,
    logistic_transform,
    normalize_sampling_weights,
    validate_single_point_scores_df,
)
from local_pipeline.common.sequence_utils import mutation_human_readable  # noqa: E402

def parse_args():
    p = argparse.ArgumentParser(description='Step 02: compute mutation sampling weights')
    p.add_argument('--input', required=True, help='单点打分 CSV（见 README）')
    p.add_argument('--output', required=True, help='输出 sampling weights CSV')
    p.add_argument('--no-normalize', action='store_true',
                   help='不归一化权重（仅输出未归一化 logistic 之和）')
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    validate_single_point_scores_df(df)

    df = df.copy()
    df['location'] = df['location'].astype(int)
    df['mutationHumanReadable'] = df.apply(
        lambda r: mutation_human_readable(r['original_aa'], r['location'], r['mutant_aa']),
        axis=1,
    )
    for col in SCORE_COLUMNS:
        df['l_' + col] = df[col].apply(
            lambda v: logistic_transform(v) if pd.notna(v) else 0.0
        )
    df['samplingWeight'] = df.apply(combine_tool_scores, axis=1)

    if not args.no_normalize:
        df = normalize_sampling_weights(df)

    out_cols = [
        'mutation', 'mutationHumanReadable', 'chain', 'location',
        'original_aa', 'mutant_aa', 'samplingWeight',
    ] + SCORE_COLUMNS + ['l_' + c for c in SCORE_COLUMNS]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df[out_cols].to_csv(args.output, index=False)
    print('Wrote {} rows to {}'.format(len(df), args.output))


if __name__ == '__main__':
    main()
