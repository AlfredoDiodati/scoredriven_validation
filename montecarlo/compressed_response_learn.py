# Learns low-dimensional compressions of every simulated run's impulse
# responses, for montecarlo/compressed_response_recovery.c to score.
#
# Reads the response caches montecarlo/sweep_grid.c keeps in
# out/sweep_grid_response_cache/ (one float32 vector per run: 525 entries for
# the t-QVARMA, 400 for the linear local projection). Everything is learned
# from runs 500-999 of every configuration; runs 0-499 are never touched, so
# a benchmark drawn from them is new data to the compression.
#
# Two compressions per auxiliary model, and per stack of models:
#   discriminant_<k>: the k directions along which the configurations' average
#     responses spread most relative to the spread of one configuration's runs
#     around its own average (generalised eigenvectors of the between- against
#     the pooled within-configuration covariance), scaled so the within-
#     configuration variance is one in every direction.
#   parameters: ridge regression of the nine standardised design parameters
#     on the response vector; the compression is the nine fitted values.
#
# Writes the compressed vectors of all runs, in the response-cache layout, to
# out/compressed_response_cache/<label>.f32 (ignored by git), the raw caches
# restricted to runs 0-199 beside them, and what was learned to
# montecarlo/out/compressed_response/learn_<label>.txt, <label> being the model
# names joined by underscores.
#
# usage: python montecarlo/compressed_response_learn.py qvarma
#        python montecarlo/compressed_response_learn.py qvarma+lp_lin
#        python montecarlo/compressed_response_learn.py lp_nl raw_only
# raw_only writes the raw cache restricted to runs 0-199 and learns nothing.
# DIRECTIONS (default 3,5,8,12,20,40) and RIDGE_RUNS (training runs per
# configuration in the regression, default 40) override the defaults.
import os
import sys

import numpy as np
import pandas as pd
import scipy.linalg as sl

RESPONSE_CACHE = 'out/sweep_grid_response_cache'
CACHE_DIR = 'out/compressed_response_cache'
REPORT_DIR = 'montecarlo/out/compressed_response'
TRAINING_FROM = 500
RAW_RUNS_KEPT = 200
RIDGE_RUNS = int(os.environ.get('RIDGE_RUNS', '40'))
DIRECTIONS = [int(k) for k in os.environ.get('DIRECTIONS', '3,5,8,12,20,40').split(',')]


def open_cache(name):
    _, n_cops, n_runs, dim = np.fromfile(f'{RESPONSE_CACHE}/{name}.f32', dtype=np.int32, count=4)
    return np.memmap(f'{RESPONSE_CACHE}/{name}.f32', dtype=np.float32, mode='r', offset=16,
                     shape=(n_cops, n_runs, dim))


def write_cache(path, array):
    # Written beside and renamed, so a reader holding the old file keeps it.
    with open(path + '.partial', 'wb') as f:
        np.array([99, *array.shape], dtype=np.int32).tofile(f)
        np.ascontiguousarray(array, dtype=np.float32).tofile(f)
    os.replace(path + '.partial', path)


def write_raw_subset(name):
    """The raw cache restricted to runs 0-199, streamed per configuration,
    so the baseline can be scored without holding 2 GB in memory."""
    cache = open_cache(name)
    n_cops, _, dim = cache.shape
    path = f'{CACHE_DIR}/{name}_raw.f32'
    with open(path + '.partial', 'wb') as f:
        np.array([99, n_cops, RAW_RUNS_KEPT, dim], dtype=np.int32).tofile(f)
        for c in range(n_cops):
            np.asarray(cache[c, :RAW_RUNS_KEPT]).tofile(f)
    os.replace(path + '.partial', path)


def training_moments(caches):
    """Each configuration's average response over the training runs and the
    pooled within-configuration covariance, for the stacked models."""
    n_cops = caches[0].shape[0]
    dim = sum(cache.shape[2] for cache in caches)
    mean = np.zeros((n_cops, dim))
    within = np.zeros((dim, dim))
    degrees = 0
    for c in range(n_cops):
        X = np.hstack([np.asarray(cache[c, TRAINING_FROM:], dtype=np.float64) for cache in caches])
        mean[c] = X.mean(0)
        Z = X - mean[c]
        within += Z.T @ Z
        degrees += len(Z) - 1
    return mean, within / degrees


def discriminant_directions(mean, within):
    # Entries no run moves, such as an impact response fixed at one by the
    # unit shock, carry no information and would make within singular.
    live = np.diag(within) > 1e-12 * np.diag(within).max()
    W = within[np.ix_(live, live)]
    W = W + 1e-6 * np.trace(W) / len(W) * np.eye(len(W))
    B = np.cov(mean[:, live].T)
    ratio, V = sl.eigh(B, W)
    order = np.argsort(ratio)[::-1]
    return live, V[:, order], ratio[order]


def project_all(caches, transform, dim_out):
    n_cops, n_runs = caches[0].shape[:2]
    out = np.zeros((n_cops, n_runs, dim_out), dtype=np.float32)
    for c in range(n_cops):
        X = np.hstack([np.asarray(cache[c], dtype=np.float64) for cache in caches])
        out[c] = transform(X)
    return out


def parameter_regression(caches, live, theta, report):
    """Ridge regression of the standardised parameters on the training runs,
    penalty chosen on a held-out tenth of the configurations."""
    rng = np.random.default_rng(0)
    n_cops = caches[0].shape[0]
    rows = []
    for c in range(n_cops):
        pick = np.sort(TRAINING_FROM + rng.choice(caches[0].shape[1] - TRAINING_FROM, RIDGE_RUNS, replace=False))
        rows.append(np.hstack([np.asarray(cache[c, pick], dtype=np.float64) for cache in caches])[:, live])
    X = np.vstack(rows)
    del rows
    Y = np.repeat(theta, RIDGE_RUNS, axis=0)
    centre, scale = X.mean(0), X.std(0) + 1e-12
    X -= centre
    X /= scale
    held_out = np.repeat(rng.permutation(n_cops) < n_cops // 10, RIDGE_RUNS)
    best_error, best_penalty = np.inf, None
    for penalty in (1e-1, 1e0, 1e1, 1e2, 1e3, 1e4):
        A = X[~held_out]
        coefficients = sl.solve(A.T @ A + penalty * np.eye(A.shape[1]), A.T @ Y[~held_out], assume_a='pos')
        error = ((X[held_out] @ coefficients - Y[held_out]) ** 2).mean(0)
        report.write(f'ridge penalty {penalty:g}: held-out R2 per parameter {np.round(1 - error, 2).tolist()}\n')
        if error.mean() < best_error:
            best_error, best_penalty = error.mean(), penalty
    report.write(f'chosen penalty {best_penalty:g}\n')
    coefficients = sl.solve(X.T @ X + best_penalty * np.eye(X.shape[1]), X.T @ Y, assume_a='pos')
    return centre, scale, coefficients


def main():
    names = sys.argv[1].split('+')
    label = '_'.join(names)
    os.makedirs(CACHE_DIR, exist_ok=True)
    if sys.argv[2:] == ['raw_only']:
        write_raw_subset(label)
        return
    os.makedirs(REPORT_DIR, exist_ok=True)
    design = pd.read_csv('dataset/abm_system_design.csv')
    theta = ((design - design.mean()) / design.std(ddof=0)).to_numpy()
    caches = [open_cache(name) for name in names]

    with open(f'{REPORT_DIR}/learn_{label}.txt', 'w') as report:
        report.write(f'{label}: training runs {TRAINING_FROM}-{caches[0].shape[1] - 1} of every configuration\n')
        if len(names) == 1 and not os.path.exists(f'{CACHE_DIR}/{label}_raw.f32'):
            write_raw_subset(names[0])

        mean, within = training_moments(caches)
        live, V, ratio = discriminant_directions(mean, within)
        report.write(f'entries {len(live)}, of which some run moves {live.sum()}\n')
        report.write('between- over within-configuration variance, leading 30 directions: '
                     f'{np.round(ratio[:30], 2).tolist()}\n')
        report.write(f'directions with ratio above 1: {(ratio > 1).sum()}\n')
        for k in DIRECTIONS:
            Vk = V[:, :k]
            write_cache(f'{CACHE_DIR}/{label}_discriminant_{k}.f32', project_all(caches, lambda X: X[:, live] @ Vk, k))

        centre, scale, coefficients = parameter_regression(caches, live, theta, report)
        write_cache(f'{CACHE_DIR}/{label}_parameters.f32',
                    project_all(caches, lambda X: ((X[:, live] - centre) / scale) @ coefficients, theta.shape[1]))


if __name__ == '__main__':
    main()
