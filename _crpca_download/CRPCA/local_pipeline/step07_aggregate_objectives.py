#!/usr/bin/env python3
"""Step 07 — 聚合多点 objective。"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.objective_aggregate import aggregate_objectives  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 07: aggregate multipoint objectives')
    p.add_argument('--sequences', required=True)
    p.add_argument('--single-point-scores', required=True)
    p.add_argument('--output', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    seq_df = pd.read_csv(args.sequences)
    sp = pd.read_csv(args.single_point_scores)
    out = aggregate_objectives(seq_df, sp)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print('Wrote {} rows to {}'.format(len(out), args.output))


if __name__ == '__main__':
    main()
