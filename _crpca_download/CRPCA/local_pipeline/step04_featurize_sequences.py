#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 04 — 对序列计算 86 维界面特征（化学/大小类别配对计数 + WT 对照）。

用法:
  python step04_featurize_sequences.py \\
    --sequences sequences/generated.csv \\
    --master-fasta configs/master.fasta \\
    --interface-pairs structure/interface_pairs.json \\
    --output features/generated_features.csv
"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.featurization import (  # noqa: E402
    featurize_sequences_df,
    load_interface_pairs,
)


def parse_args():
    p = argparse.ArgumentParser(description='Step 04: 86-d interface featurization')
    p.add_argument('--sequences', required=True, help='Step 03 输出 CSV（含 sequence 列）')
    p.add_argument('--master-fasta', required=True)
    p.add_argument('--interface-pairs', required=True, help='Step 03b JSON')
    p.add_argument('--output', required=True)
    return p.parse_args()


def read_master_fasta(path):
    from Bio import SeqIO
    return str(next(SeqIO.parse(path, 'fasta')).seq)


def main():
    args = parse_args()
    seq_df = pd.read_csv(args.sequences)
    if 'sequence' not in seq_df.columns:
        raise ValueError('输入需包含 sequence 列')

    master = read_master_fasta(args.master_fasta)
    pairs, meta = load_interface_pairs(args.interface_pairs)
    chain = meta.get('chain_designator', '')

    feat_df = featurize_sequences_df(seq_df, master, pairs, chain_designator=chain)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    feat_df.to_csv(args.output, index=False)
    n_feat = len([c for c in feat_df.columns if c not in (
        'sequence', 'num_mutations', 'mutation_str')])
    print('Wrote {} rows, {} feature columns to {}'.format(len(feat_df), n_feat, args.output))


if __name__ == '__main__':
    main()
