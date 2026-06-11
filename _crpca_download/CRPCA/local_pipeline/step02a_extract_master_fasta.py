#!/usr/bin/env python3
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
Step 03a — 从 PDB 抗体链自动提取 master FASTA（供 --master-fasta 使用）。

将指定抗体链按顺序拼接为一条序列（默认重链在前、轻链在后），
并输出 PDB 残基编号 → 线性 1-based 位置的映射表，便于核对 CSV 中的 location。

用法:
  python step03a_extract_master_fasta.py \\
    --pdb structure/7l7e.pdb \\
    --antibody-chains M,N \\
    --output configs/master.fasta \\
    --mapping-output configs/master_residue_mapping.csv
"""

from __future__ import division, print_function

import argparse
import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(_ROOT.parent))


def parse_args():
    p = argparse.ArgumentParser(description='Step 03a: extract master FASTA from PDB')
    p.add_argument('--pdb', required=True, help='PDB/mmCIF 文件路径')
    p.add_argument('--antibody-chains', required=True,
                   help='逗号分隔；按顺序拼接，如 M,N 表示重链+轻链')
    p.add_argument('--fasta-id', default='master_from_pdb',
                   help='FASTA header 名称（不影响下游计算）')
    p.add_argument('--output', required=True, help='输出 master.fasta')
    p.add_argument('--mapping-output', default=None,
                   help='可选：残基编号映射 CSV')
    p.add_argument('--include-hetero', action='store_true',
                   help='包含 HETATM 标准氨基酸（默认仅标准残基）')
    return p.parse_args()


def pdb_residue_id(residue):
    """PDB 残基编号字符串（含 insertion code）。"""
    het, resseq, icode = residue.id
    num = str(resseq)
    if icode and str(icode).strip():
        num += str(icode).strip()
    return num


def chain_to_sequence(chain, include_hetero=False):
    """从一条链提取氨基酸序列及残基元数据。"""
    seq_chars = []
    rows = []
    for residue in chain:
        het, resseq, icode = residue.id
        if het != ' ' and not include_hetero:
            continue
        if not is_aa(residue, standard=True):
            continue
        if 'CA' not in residue:
            continue
        try:
            aa = three_to_one(residue.get_resname())
        except KeyError:
            continue
        seq_chars.append(aa)
        rows.append({
            'chain': chain.id,
            'pdb_resnum': pdb_residue_id(residue),
            'resseq': resseq,
            'icode': str(icode).strip(),
            'aa': aa,
            'linear_index': len(seq_chars),  # 1-based position in master
        })
    return ''.join(seq_chars), rows


def write_fasta(path, seq_id, sequence):
    with open(path, 'w') as f:
        f.write('>{}\n'.format(seq_id))
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + '\n')


def main():
    args = parse_args()
    chain_ids = [c.strip() for c in args.antibody_chains.split(',') if c.strip()]
    if not chain_ids:
        raise ValueError('--antibody-chains 不能为空')

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('cpx', args.pdb)
    model = structure[0]
    pdb_chains = {chain.id: chain for chain in model}

    missing = [c for c in chain_ids if c not in pdb_chains]
    if missing:
        raise ValueError('PDB 中找不到链: {}；现有链: {}'.format(
            missing, sorted(pdb_chains.keys())))

    parts = []
    mapping_rows = []
    offset = 0
    for cid in chain_ids:
        seq_part, rows = chain_to_sequence(pdb_chains[cid], args.include_hetero)
        if not seq_part:
            raise ValueError('链 {} 未提取到氨基酸序列'.format(cid))
        for r in rows:
            r['linear_index'] = r['linear_index'] + offset
        offset += len(seq_part)
        parts.append(seq_part)
        mapping_rows.extend(rows)

    master_seq = ''.join(parts)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_path, args.fasta_id, master_seq)

    print('Wrote master FASTA: {} ({} aa, chains: {})'.format(
        out_path, len(master_seq), ','.join(chain_ids)))

    if args.mapping_output:
        import pandas as pd
        map_path = Path(args.mapping_output)
        map_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(mapping_rows).to_csv(map_path, index=False)
        print('Wrote residue mapping: {} ({} residues)'.format(
            map_path, len(mapping_rows)))


if __name__ == '__main__':
    main()
