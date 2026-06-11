#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Step 06 — MEI 选序列（GP 不重训，供 BO 循环调用）。"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.mei_selection import mei_select_sequences  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 06: MEI sequence selection')
    p.add_argument('--features', required=True)
    p.add_argument('--model-dir', required=True)
    p.add_argument('--single-point-scores', default=None,
                   help='用于 best_so_far 的 rosetta_flex 最小值')
    p.add_argument('--batch-size', type=int, default=50)
    p.add_argument('--output', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.features)
    scores = pd.read_csv(args.single_point_scores) if args.single_point_scores else None
    out, best_sf = mei_select_sequences(df, args.model_dir, args.batch_size, scores)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print('Wrote top {} sequences by MEI to {}'.format(len(out), args.output))
    print('best_so_far (most negative ddG): {:.4f}'.format(best_sf))


if __name__ == '__main__':
    main()
