# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""MLP + single-task GP for Rosetta Flex ddG (paper Step 2.5)."""

from __future__ import division, print_function

import gpytorch
import torch


class DKLExactGP(gpytorch.models.ExactGP):
    """86-dim features -> Linear(86,40) -> tanh -> Linear(40,10) -> RBF GP."""

    def __init__(self, train_x, train_y, likelihood, mlp_hidden=40, mlp_out=10):
        super(DKLExactGP, self).__init__(train_x, train_y, likelihood)
        self.feature_extractor = torch.nn.Sequential(
            torch.nn.Linear(train_x.size(-1), mlp_hidden),
            torch.nn.Tanh(),
            torch.nn.Linear(mlp_hidden, mlp_out),
        )
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())

    def forward(self, x):
        projected = self.feature_extractor(x)
        mean_x = self.mean_module(projected)
        covar_x = self.covar_module(projected)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def train_gp(
    x_tensor,
    y_tensor,
    num_iters=500,
    lr=0.01,
    mlp_hidden=40,
    mlp_out=10,
    feature_scale=1.0,
):
    """
    Fit GP on (scaled) features and Rosetta ddG targets.

    Returns trained model, likelihood, training loss history.
    """
    x = x_tensor.float() * feature_scale
    y = y_tensor.float()
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = DKLExactGP(x, y, likelihood, mlp_hidden=mlp_hidden, mlp_out=mlp_out)
    model.train()
    likelihood.train()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
    losses = []
    for i in range(num_iters):
        optimizer.zero_grad()
        output = model(x)
        loss = -mll(output, y)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())
        if (i + 1) % 100 == 0:
            print('Iter {}/{}  -mll={:.4f}'.format(i + 1, num_iters, loss.item()))
    return model, likelihood, losses


def predict_gp(model, likelihood, x_tensor, feature_scale=1.0):
    model.eval()
    likelihood.eval()
    x = x_tensor.float() * feature_scale
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        pred = likelihood(model(x))
    return pred.mean, pred.stddev


def compute_mei_scores(mean, stddev, best_so_far):
    """MEI for minimization (more negative ddG is better). best_so_far is most negative seen."""
    a = torch.distributions.Normal(0.0, 1.0)
    z = (best_so_far - mean) / stddev.clamp(min=1e-6)
    z = torch.clamp(z, -6.0, 6.0)
    scores = (best_so_far - mean) * a.cdf(z) + stddev * torch.exp(a.log_prob(z))
    return scores
