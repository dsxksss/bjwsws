#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 03b（可选）— 从 PDB 提取 WT 界面残基对，供 Step 04 特征化使用。

判定：抗体链 Cα 与抗原链 Cα 距离 < cutoff（默认 10 Å）。

用法:
  python step03b_extract_interface_pairs.py \\
    --pdb structure/7l7e.pdb \\
    --antibody-chains M,N \\
    --antigen-chains S \\
    --cutoff 10.0 \\
    --output structure/interface_pairs.json
"""

from __future__ import division, print_function

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))

from local_pipeline.common.pdb_structure import resolve_antigen_chains  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Step 03b: extract interface pairs from PDB')
    p.add_argument('--pdb', required=True)
    p.add_argument('--antibody-chains', required=True, help='逗号分隔，可变链（抗体）')
    p.add_argument('--antigen-chains', default=None,
                   help='逗号分隔，上下文链（抗原）；省略则取 PDB 中除抗体链外的所有链')
    p.add_argument('--cutoff', type=float, default=10.0)
    p.add_argument('--master-seq', default=None, help='可选，写入 JSON 供下游校验')
    p.add_argument('--output', required=True)
    return p.parse_args()


def ca_coord(residue):
    return residue['CA'].get_coord()


def residue_triple(residue):
    rid = residue.full_id[3]
    num = ''.join(str(x) for x in rid if str(x).strip())
    return (residue.full_id[2], num, three_to_one(residue.resname))


def extract_pairs(pdb_path, ab_chains, ag_chains, cutoff):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('cpx', pdb_path)
    model = structure[0]
    ab_res = []
    ag_res = []
    for chain in model:
        if chain.id in ab_chains:
            for res in chain:
                if 'CA' in res:
                    ab_res.append(res)
        elif chain.id in ag_chains:
            for res in chain:
                if 'CA' in res:
                    ag_res.append(res)
    pairs = []
    for ra in ab_res:
        for rb in ag_res:
            d = np.linalg.norm(ca_coord(ra) - ca_coord(rb))
            if d < cutoff:
                pairs.append((residue_triple(ra), residue_triple(rb)))
    pairs = sorted(set(pairs))
    return pairs


def main():
    args = parse_args()
    ab = [c.strip() for c in args.antibody_chains.split(',') if c.strip()]
    ag = resolve_antigen_chains(args.pdb, ab, args.antigen_chains)
    pairs = extract_pairs(args.pdb, ab, ag, args.cutoff)
    pair_list = [[list(a), list(b)] for a, b in pairs]
    payload = {
        'pdb': str(Path(args.pdb).resolve()),
        'antibody_chains': ab,
        'antigen_chains': ag,
        'cutoff_angstrom': args.cutoff,
        'pairs': pair_list,
        'chain_designator': '',
    }
    if args.master_seq:
        payload['master_seq'] = args.master_seq
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(payload, f, indent=2)
    print('Wrote {} interface pairs to {}'.format(len(pair_list), args.output))


if __name__ == '__main__':
    main()
