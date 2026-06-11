# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""Library functions for MEI selection (Step 06)."""

from __future__ import division, print_function

import json
import pickle
from pathlib import Path

import pandas as pd
import torch
import gpytorch

from local_pipeline.common.gp_model import DKLExactGP, compute_mei_scores, predict_gp


def load_gp_model(model_dir):
    model_dir = Path(model_dir)
    with open(model_dir / 'model_meta.json', 'r') as f:
        meta = json.load(f)
    with open(model_dir / 'training_data.pkl', 'rb') as f:
        train = pickle.load(f)
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = DKLExactGP(
        train['x'], train['y'], likelihood,
        mlp_hidden=meta['mlp_hidden'],
        mlp_out=meta['mlp_out'],
    )
    model.load_state_dict(torch.load(model_dir / 'gp_model_state.pth', map_location='cpu'))
    likelihood.load_state_dict(
        torch.load(model_dir / 'gp_likelihood_state.pth', map_location='cpu'))
    return model, likelihood, meta


def mei_select_sequences(features_df, model_dir, batch_size, single_point_scores_df=None):
    """
    :return: DataFrame top batch_size by MEI, sorted descending
    """
    model, likelihood, meta = load_gp_model(model_dir)
    pred_cols = meta['predictor_columns']
    missing = [c for c in pred_cols if c not in features_df.columns]
    if missing:
        raise ValueError('特征表缺少列: {}'.format(missing))

    x = torch.tensor(features_df[pred_cols].values, dtype=torch.float)
    mean, std = predict_gp(model, likelihood, x, feature_scale=meta['feature_scale'])

    if single_point_scores_df is not None and 'rosetta_flex' in single_point_scores_df.columns:
        best_so_far = float(single_point_scores_df['rosetta_flex'].min())
    else:
        with open(Path(model_dir) / 'training_data.pkl', 'rb') as f:
            best_so_far = float(pickle.load(f)['y'].min().item())

    mei = compute_mei_scores(mean, std, best_so_far)
    out = features_df.copy()
    out['pred_mean'] = mean.numpy()
    out['pred_std'] = std.numpy()
    out['mei_score'] = mei.numpy()
    out = out.sort_values('mei_score', ascending=False).head(batch_size)
    return out, best_so_far
