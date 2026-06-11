# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Shared I/O helpers for the local pipeline."""

from __future__ import division, print_function

from Bio import SeqIO


def read_master_fasta(path):
    return str(next(SeqIO.parse(path, 'fasta')).seq)


def write_fasta(path, seq_id, sequence):
    with open(path, 'w') as f:
        f.write('>{}\n'.format(seq_id))
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + '\n')


def write_recommendations_fasta(df, path, id_prefix='candidate'):
    with open(path, 'w') as f:
        for i, row in df.iterrows():
            muts = row.get('mutationHumanReadable', 'unknown')
            f.write('>{}_{} {}\n'.format(id_prefix, i, muts))
            seq = row['sequence']
            for j in range(0, len(seq), 80):
                f.write(seq[j:j + 80] + '\n')
