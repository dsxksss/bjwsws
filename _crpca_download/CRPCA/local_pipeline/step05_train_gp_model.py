#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Step 05 — 训练 GP（仅执行一次，BO 循环中不重训）。"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.gp_training import train_gp_model  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 05: train MLP+GP surrogate (once)')
    p.add_argument('--single-point-features', required=True)
    p.add_argument('--single-point-scores', required=True)
    p.add_argument('--target-column', default='rosetta_flex')
    p.add_argument('--join-on', default='mutationHumanReadable')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--num-iters', type=int, default=500)
    p.add_argument('--lr', type=float, default=0.01)
    p.add_argument('--feature-scale', type=float, default=1.0)
    p.add_argument('--mlp-hidden', type=int, default=40)
    p.add_argument('--mlp-out', type=int, default=10)
    return p.parse_args()


def main():
    args = parse_args()
    feat_df = pd.read_csv(args.single_point_features)
    scores_df = pd.read_csv(args.single_point_scores)
    meta = train_gp_model(
        feat_df, scores_df, args.output_dir,
        target_column=args.target_column,
        join_on=args.join_on,
        num_iters=args.num_iters,
        lr=args.lr,
        feature_scale=args.feature_scale,
        mlp_hidden=args.mlp_hidden,
        mlp_out=args.mlp_out,
    )
    print('Saved model to {} (n_train={})'.format(args.output_dir, meta['n_train']))


if __name__ == '__main__':
    main()
