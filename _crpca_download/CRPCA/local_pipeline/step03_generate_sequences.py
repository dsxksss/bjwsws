#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 03 — 序列生成器（论文 2.3 完整版：Step A + Step B）。

Step A: k ~ Uniform[min_locations, max_locations]
Step B: 按 samplingWeight 无放回加权采样 k 个突变
"""

from __future__ import division, print_function

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.mutation_generator import generate_mutant_sequences  # noqa: E402
from local_pipeline.common.pipeline_io import read_master_fasta  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 03: sequence generator (Step A + B)')
    p.add_argument('--master-fasta', required=True)
    p.add_argument('--allowed-mutations', required=True)
    p.add_argument('--sampling-weights', required=True)
    p.add_argument('--min-locations', type=int, default=1, help='Step A 最少突变位点')
    p.add_argument('--max-locations', type=int, default=8, help='Step A 最多突变位点')
    p.add_argument('--num-sequences', type=int, default=1000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--output', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    master = read_master_fasta(args.master_fasta)
    with open(args.allowed_mutations, 'r') as f:
        allowed = json.load(f)
    weights = pd.read_csv(args.sampling_weights)

    rows = generate_mutant_sequences(
        master, weights, allowed,
        number_to_generate=args.num_sequences,
        min_locations=args.min_locations,
        max_locations=args.max_locations,
        seed=args.seed,
    )
    out = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print('Wrote {} sequences to {}'.format(len(out), args.output))


if __name__ == '__main__':
    main()
