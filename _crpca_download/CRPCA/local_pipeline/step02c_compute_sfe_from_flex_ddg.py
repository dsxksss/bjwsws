#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 02c — 从 21×2 个 Rosetta Flex ddG 聚合 SFE 分数。

论文公式:
  SFE = (mean(forward, IQR-filtered) - mean(reverse, IQR-filtered)) / 2

输入宽表 CSV — examples/sfe_flex_ddg.example.csv:
  mutation,chain,forward_01,...,forward_21,reverse_01,...,reverse_21

用法:
  python step02c_compute_sfe_from_flex_ddg.py \\
    --input examples/sfe_flex_ddg.example.csv \\
    --output work/sfe_scores.csv
"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.sfe_aggregation import (  # noqa: E402
    DEFAULT_N_CONFORMATIONS,
    aggregate_sfe_table,
    load_flex_ddg_wide,
)


def parse_args():
    p = argparse.ArgumentParser(
        description='Aggregate SFE ddG from Rosetta Flex forward/reverse ensembles')
    p.add_argument('--input', required=True, help='宽表 CSV')
    p.add_argument('--output', required=True, help='输出 CSV')
    p.add_argument(
        '--n-conformations', type=int, default=DEFAULT_N_CONFORMATIONS,
        help='每组构象数（默认 21）')
    p.add_argument(
        '--iqr-multiplier', type=float, default=1.5,
        help='Tukey outlier 阈值倍数（默认 1.5）')
    p.add_argument(
        '--mode', choices=['tukey', 'strict-iqr'], default='tukey',
        help='outlier 过滤模式')
    p.add_argument(
        '--invert-for-sampling', action='store_true',
        help='输出 sfe = -sfe_ddg（越大越 desirable，供 2.2 CSV 使用）')
    return p.parse_args()


def main():
    args = parse_args()

    records = load_flex_ddg_wide(
        args.input, n_conformations=args.n_conformations)

    out = aggregate_sfe_table(
        records,
        iqr_multiplier=args.iqr_multiplier,
        mode=args.mode,
        n_conformations=args.n_conformations,
        invert_for_sampling=args.invert_for_sampling,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    print('Wrote {} mutations to {}'.format(len(out), out_path))
    print('conformations: {} | mode: {}'.format(args.n_conformations, args.mode))
    if len(out):
        r0 = out.iloc[0]
        print('Example: {} -> sfe_ddg={:.4f}, forward_mean={:.4f}, '
              'reverse_mean={:.4f}'.format(
                  r0['mutation'], r0['sfe_ddg'],
                  r0['forward_mean'], r0['reverse_mean']))


if __name__ == '__main__':
    main()
