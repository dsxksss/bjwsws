# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
论文 2.3 序列生成器：Step A + Step B（与 mutant_generator_sampling.py 一致）。

Step A: 突变位点数 k ~ Uniform[min_locations, max_locations]
Step B: 按 samplingWeight 无放回加权采样 k 个单点突变
"""

from __future__ import division, print_function

import random

import numpy as np
import pandas as pd

from local_pipeline.common.sequence_utils import mutate_seq


def prepare_sampling_table(weights_df, allowed_hr):
    df = weights_df.copy()
    if 'mutationHumanReadable' not in df.columns:
        df['mutationHumanReadable'] = df.apply(
            lambda r: '{}{}{}'.format(r['original_aa'], int(r['location']), r['mutant_aa']),
            axis=1,
        )
    df['mutation'] = df.apply(
        lambda r: ('', str(int(r['location'])), r['original_aa'], r['mutant_aa']),
        axis=1,
    )
    df['location'] = df['location'].astype(str)
    df = df.loc[df['mutationHumanReadable'].isin(allowed_hr)]
    if df.empty:
        raise ValueError('采样表与 allowed_mutations 无交集')
    return df


def _sample_mutations_for_one(data, num_locations, max_tries=200):
    mutations = set()
    tries = 0
    while len(mutations) < num_locations:
        tries += 1
        if tries > max_tries:
            raise ValueError(
                '无法在 {} 次尝试内采样 {} 个不重复位点'.format(max_tries, num_locations)
            )
        pick = random.choices(
            list(data['mutation']), weights=list(data['samplingWeight']), k=1
        )[0]
        if pick[2] == pick[3]:
            continue
        if pick[1] in [m[1] for m in mutations]:
            continue
        mutations.add(pick)
    return list(mutations)


def generate_mutant_sequences(
        master_sequence,
        sampling_weights_df,
        allowed_mutations,
        number_to_generate,
        min_locations=1,
        max_locations=8,
        exclude_sequences=None,
        max_tries_per_mutant=200,
        seed=None,
):
    """
    :param allowed_mutations: list of [location, original_aa, [to_aa,...]]
    :param exclude_sequences: set of sequences to skip (跨 BO 轮去重)
    :return: list of dict with sequence, num_mutations, mutations, mutationHumanReadable
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    allowed_hr = set()
    for loc, orig, to_list in allowed_mutations:
        for to_aa in to_list:
            allowed_hr.add('{}{}{}'.format(orig, loc, to_aa))

    data = prepare_sampling_table(sampling_weights_df, allowed_hr)
    locations = set(data['location'])
    if len(locations) < min_locations:
        raise ValueError(
            '可用位点数 {} < min_locations {}'.format(len(locations), min_locations)
        )
    max_locations = min(max_locations, len(locations))

    if exclude_sequences is None:
        exclude_sequences = set()

    rows = []
    seen = exclude_sequences.copy()
    attempts = 0
    max_attempts = max(number_to_generate * 20, number_to_generate + 1)
    while len(rows) < number_to_generate and attempts < max_attempts:
        attempts += 1
        k = int(np.random.randint(min_locations, max_locations + 1))
        muts = _sample_mutations_for_one(data, k, max_tries=max_tries_per_mutant)
        seq = mutate_seq(master_sequence, muts)
        if seq in seen:
            continue
        seen.add(seq)
        rows.append({
            'sequence': seq,
            'num_mutations': len(muts),
            'mutations': str(muts),
            'mutationHumanReadable': ','.join(
                ['{}{}{}'.format(m[2], m[1], m[3]) for m in muts]
            ),
        })
    return rows
