#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 08 — Pareto 非支配筛选 + 加权降维（简化版，无 HPC）。

用法:
  python step08_pareto_select.py \\
    --objectives scores/multipoint_objectives.csv \\
    --objective-columns sum_rosetta_flex,sum_foldx,sum_abnativ,num_mutations \\
    --maximize-columns sum_abnativ \\
    --downselect-quota 100 \\
    --output selection/pareto_final.csv
"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

import abag_ml.pareto_selection as ps  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 08: Pareto + weighted downselect')
    p.add_argument('--objectives', required=True)
    p.add_argument('--objective-columns', required=True,
                   help='逗号分隔；默认越小越好')
    p.add_argument('--maximize-columns', default='',
                   help='逗号分隔，需最大化的列（如 abnativ）')
    p.add_argument('--epsilons', default='',
                   help='逗号分隔，与 objective-columns 对齐')
    p.add_argument('--downselect-quota', type=int, default=0,
                   help='>0 时从 Pareto 集做加权降维')
    p.add_argument('--weights', default='',
                   help='降维权重，逗号分隔')
    p.add_argument('--output', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.objectives)
    obj_cols = [c.strip() for c in args.objective_columns.split(',')]
    max_cols = set(c.strip() for c in args.maximize_columns.split(',') if c.strip())
    eps = [float(x) for x in args.epsilons.split(',')] if args.epsilons else [0.0] * len(obj_cols)
    if len(eps) != len(obj_cols):
        raise ValueError('epsilons 数量需与 objective-columns 一致')

    dom_funcs = []
    for col in obj_cols:
        if col in max_cols:
            dom_funcs.append(lambda r=None, c=col: ps.simple_negative_row_scorer(r, column=c))
        else:
            dom_funcs.append(lambda r=None, c=col: ps.simple_row_scorer(r, column=c))

    pareto_idx = ps.get_pareto_rows(
        df, dominance_functions=dom_funcs, scalar_epsilons=eps, returnints=True
    )
    df['ParetoSet'] = [i in pareto_idx for i in range(len(df))]

    if args.downselect_quota > 0:
        if not args.weights:
            weights = [-1.0] * len(obj_cols)
        else:
            weights = [float(w) for w in args.weights.split(',')]
        sel_cols = obj_cols
        df['downselection_score'] = 0.0
        for col, w in zip(sel_cols, weights):
            if col in max_cols:
                df['downselection_score'] += w * df[col]
            else:
                df['downselection_score'] += w * df[col]
        tmp = df['downselection_score'].copy()
        tmp[~df['ParetoSet']] = np.nan
        best = []
        for _ in range(args.downselect_quota):
            try:
                i = np.nanargmax(tmp.values)
            except ValueError:
                break
            best.append(i)
            tmp.iloc[i] = np.nan
        df['DownSelected'] = [i in best for i in range(len(df))]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print('Pareto set size: {}'.format(sum(df['ParetoSet'])))
    if args.downselect_quota > 0:
        print('DownSelected: {}'.format(sum(df.get('DownSelected', []))))


if __name__ == '__main__':
    main()
