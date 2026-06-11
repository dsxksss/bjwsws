#!/usr/bin/env python3
"""从单点打分表生成单点全长序列 CSV，供 Step 04/05 训练 GP。"""

from __future__ import division, print_function

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.sequence_utils import mutate_seq  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--master-fasta', required=True)
    p.add_argument('--single-point-scores', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()

    from Bio import SeqIO
    master = str(next(SeqIO.parse(args.master_fasta, 'fasta')).seq)
    df = pd.read_csv(args.single_point_scores)
    rows = []
    for _, r in df.iterrows():
        mut = ('', str(int(r['location'])), r['original_aa'], r['mutant_aa'])
        seq = mutate_seq(master, [mut])
        hr = '{}{}{}'.format(r['original_aa'], int(r['location']), r['mutant_aa'])
        rows.append({'sequence': seq, 'mutationHumanReadable': hr, 'num_mutations': 1})
    out = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print('Wrote {} single-point sequences'.format(len(out)))


if __name__ == '__main__':
    main()
