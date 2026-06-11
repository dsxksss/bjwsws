# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Library functions for GP training (Step 05)."""

from __future__ import division, print_function

import json
import pickle
from pathlib import Path

import pandas as pd
import torch

from local_pipeline.common.gp_model import train_gp

META_COLS = {
    'sequence', 'num_mutations', 'mutation_str', 'mutationHumanReadable',
    'mutation', 'bo_round',
}


def infer_predictor_columns(df):
    skip = set(META_COLS) | {
        'chain', 'location', 'original_aa', 'mutant_aa', 'samplingWeight',
        'sfe', 'fep', 'rosetta_flex', 'foldx', 'abnativ',
        'pred_mean', 'pred_std', 'mei_score',
    }
    skip |= {c for c in df.columns if c.startswith('l_')}
    cols = []
    for c in df.columns:
        if c in skip or c.startswith('wt_'):
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            cols.append(c)
    return cols


def train_gp_model(
        features_df,
        scores_df,
        output_dir,
        target_column='rosetta_flex',
        join_on='mutationHumanReadable',
        num_iters=500,
        lr=0.01,
        feature_scale=1.0,
        mlp_hidden=40,
        mlp_out=10,
):
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if join_on in features_df.columns and scores_df is not None:
        merged = features_df.merge(
            scores_df[[join_on, target_column]], on=join_on, how='inner')
    elif target_column in features_df.columns:
        merged = features_df
    else:
        raise ValueError('无法获取训练标签 {}'.format(target_column))

    pred_cols = infer_predictor_columns(merged)
    if not pred_cols:
        raise ValueError('未识别到数值特征列')

    x = torch.tensor(merged[pred_cols].values, dtype=torch.float)
    y = torch.tensor(merged[target_column].values, dtype=torch.float)

    model, likelihood, losses = train_gp(
        x, y, num_iters=num_iters, lr=lr,
        mlp_hidden=mlp_hidden, mlp_out=mlp_out,
        feature_scale=feature_scale,
    )

    torch.save(model.state_dict(), out_dir / 'gp_model_state.pth')
    torch.save(likelihood.state_dict(), out_dir / 'gp_likelihood_state.pth')
    meta = {
        'predictor_columns': pred_cols,
        'target_column': target_column,
        'feature_scale': feature_scale,
        'mlp_hidden': mlp_hidden,
        'mlp_out': mlp_out,
        'n_train': int(x.shape[0]),
        'training_loss_final': losses[-1] if losses else None,
    }
    with open(out_dir / 'model_meta.json', 'w') as f:
        json.dump(meta, f, indent=2)
    with open(out_dir / 'training_data.pkl', 'wb') as f:
        pickle.dump({'x': x, 'y': y, 'predictor_columns': pred_cols}, f)
    return meta
