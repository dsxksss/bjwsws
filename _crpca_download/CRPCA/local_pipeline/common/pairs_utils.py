# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Interface pair mutation helpers (vendored from interface_residues.py)."""

from __future__ import division, print_function


def match_residue_triples(trip, query_trip, aa_from_location=2):
    return all(
        t_i == q_i if (i != aa_from_location or q_i is not None) else True
        for i, (t_i, q_i) in enumerate(zip(trip, query_trip))
    )


def mutate_pairs_once(pairs, mutation):
    for i, pair_i in enumerate(pairs):
        match = [match_residue_triples(pair_i_j, mutation[:3]) for pair_i_j in pair_i]
        pairs[i] = tuple(
            pair_i_j[:-1] + (mutation[-1],) if m_ij else pair_i_j
            for m_ij, pair_i_j in zip(match, pair_i)
        )
    return pairs


def mutate_pairs_multiple(pairs, mutations):
    if not isinstance(mutations, list):
        mutations = list(mutations)
    for mi in mutations:
        pairs = mutate_pairs_once(pairs, mi)
    return pairs
