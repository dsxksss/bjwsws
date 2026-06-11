# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""PDB structure extraction for master FASTA and interface pairs."""

from __future__ import division, print_function

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one


def pdb_residue_id(residue):
    het, resseq, icode = residue.id
    num = str(resseq)
    if icode and str(icode).strip():
        num += str(icode).strip()
    return num


def residue_triple(residue):
    return (residue.full_id[2], pdb_residue_id(residue), three_to_one(residue.resname))


def chain_to_sequence(chain, include_hetero=False):
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
            'linear_index': len(seq_chars),
        })
    return ''.join(seq_chars), rows


def get_pdb_chain_ids(pdb_path):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure('cpx', pdb_path)[0]
    return sorted(c.id for c in model)


def infer_antigen_chains(pdb_path, antibody_chain_ids):
    ab_set = set(antibody_chain_ids)
    ag_chains = [c for c in get_pdb_chain_ids(pdb_path) if c not in ab_set]
    if not ag_chains:
        raise ValueError(
            'PDB 中除抗体链 {} 外无其他链；请显式指定 --antigen-chains'.format(
                antibody_chain_ids))
    return ag_chains


def resolve_antigen_chains(pdb_path, antibody_chain_ids, antigen_chains_arg=None):
    if antigen_chains_arg:
        return [c.strip() for c in antigen_chains_arg.split(',') if c.strip()]
    ag_chains = infer_antigen_chains(pdb_path, antibody_chain_ids)
    print('[INFO] 未指定 --antigen-chains，自动推断为: {}'.format(ag_chains))
    return ag_chains


def extract_master_from_pdb(pdb_path, antibody_chain_ids, include_hetero=False):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure('cpx', pdb_path)[0]
    pdb_chains = {c.id: c for c in model}
    missing = [c for c in antibody_chain_ids if c not in pdb_chains]
    if missing:
        raise ValueError('PDB 缺少链 {}；现有 {}'.format(missing, sorted(pdb_chains.keys())))

    parts, mapping_rows, offset = [], [], 0
    for cid in antibody_chain_ids:
        seq_part, rows = chain_to_sequence(pdb_chains[cid], include_hetero)
        if not seq_part:
            raise ValueError('链 {} 无序列'.format(cid))
        for r in rows:
            r['linear_index'] += offset
        offset += len(seq_part)
        parts.append(seq_part)
        mapping_rows.extend(rows)
    return ''.join(parts), mapping_rows


def extract_interface_pairs(pdb_path, antibody_chain_ids, antigen_chain_ids, cutoff=10.0):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure('cpx', pdb_path)[0]
    ab_res, ag_res = [], []
    for chain in model:
        if chain.id in antibody_chain_ids:
            for res in chain:
                if 'CA' in res:
                    ab_res.append(res)
        elif chain.id in antigen_chain_ids:
            for res in chain:
                if 'CA' in res:
                    ag_res.append(res)
    pairs = []
    for ra in ab_res:
        ca = ra['CA'].get_coord()
        for rb in ag_res:
            d = np.linalg.norm(ca - rb['CA'].get_coord())
            if d < cutoff:
                pairs.append((residue_triple(ra), residue_triple(rb)))
    pairs = sorted(set(pairs))
    return [[list(a), list(b)] for a, b in pairs]


def write_interface_pairs_json(pdb_path, antibody_chain_ids, antigen_chain_ids,
                               pairs, output_path, cutoff=10.0):
    payload = {
        'pdb': str(Path(pdb_path).resolve()),
        'antibody_chains': list(antibody_chain_ids),
        'antigen_chains': list(antigen_chain_ids),
        'cutoff_angstrom': cutoff,
        'pairs': pairs,
        'chain_designator': '',
    }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(payload, f, indent=2)
