# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Minimal sequence mutation helpers (standalone, no improvwf deps)."""

from __future__ import division, print_function


def diff_seqs(mutant_seq, master_seq, chain_name=''):
    """Return list of (chain, position_str, from_aa, to_aa) tuples."""
    if len(mutant_seq) != len(master_seq):
        raise ValueError('Sequences differ in length.')
    mutations = []
    for idx, (mi, wi) in enumerate(zip(mutant_seq, master_seq)):
        if mi != wi:
            mutations.append((chain_name, str(idx + 1), wi, mi))
    return mutations


def mutate_seq(master_seq, mutations, chain_name=None):
    """Apply mutation 4-tuples to master_seq; positions are 1-based in string."""
    seq = list(master_seq)
    for mut in mutations:
        chain, pos, _from, to = mut
        if chain_name is not None and chain != chain_name and chain != '':
            continue
        idx = int(pos) - 1
        seq[idx] = to
    return ''.join(seq)


def mutation_human_readable(from_aa, location, to_aa):
    return '{}{}{}'.format(from_aa, location, to_aa)


def mutations_to_tuples(mutation_str_list):
    """
    Parse mutation strings like 'G112E' or 'S32A' into 4-tuples.
    Requires original_aa in the CSV row when using structured input.
    """
    raise NotImplementedError('Use CSV columns original_aa, location, mutant_aa instead.')
