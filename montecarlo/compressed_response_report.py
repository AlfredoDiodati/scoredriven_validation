# Summarises the files montecarlo/compressed_response_recovery.c wrote into
# montecarlo/out/compressed_response/report.txt: per output vector, how often
# the benchmark's own configuration is in the confidence set and has the
# smallest mean loss, and, when it does not, how far in parameter space the
# configuration that took its place lies.
#
# File names say what each file holds: recovery_<vector>_rows<n>.csv carries
# confidence sets, rank_<vector>_rows<n>.csv only ranks; <n> is how many runs
# of each configuration were rows of the loss matrix before the benchmark's
# own run was dropped. grid_<vector>_rows1000.csv holds the rows for benchmark
# runs 0 and 1 of montecarlo/out/sweep_grid.csv.gz, the raw responses at full
# scale, in the same columns; montecarlo/compressed_response_grid_baseline.py
# copies them.
import glob
import os
import re

import numpy as np
import pandas as pd

DIR = 'montecarlo/out/compressed_response'

# The parameters the stacked responses predict with held-out R2 above one
# half (montecarlo/out/compressed_response/learn_qvarma_lp_lin.txt).
IDENTIFIED = ['psi1', 'psi3', 'alfa', 'taylor1', 'taylor2', 'taylor', 'kappa']


def neighbour_rank(design, columns):
    unit = ((design - design.min()) / (design.max() - design.min()))[columns].to_numpy()
    distance = np.sqrt(((unit[:, None] - unit[None]) ** 2).sum(-1))
    np.fill_diagonal(distance, np.inf)
    # Entry (a, b) is 0 when b is a's nearest configuration, 1 for the next.
    return np.argsort(np.argsort(distance, axis=1), axis=1)


def main():
    design = pd.read_csv('dataset/abm_system_design.csv')
    neighbour_all = neighbour_rank(design, list(design.columns))
    neighbour_identified = neighbour_rank(design, IDENTIFIED)

    rows = []
    for path in sorted(glob.glob(f'{DIR}/*_rows*.csv')):
        kind, vector, n_rows = re.match(r'(recovery|rank|grid)_(.+)_rows(\d+)\.csv', os.path.basename(path)).groups()
        d = pd.read_csv(path)
        cop, winner = d.cop.to_numpy() - 1, d.lowest_loss_cop.to_numpy() - 1
        missed = winner != cop
        row = dict(vector=vector, rows=int(n_rows), benchmarks=len(d),
                   benchmark_runs=f'{d.run.min()}-{d.run.max()}',
                   smallest_loss=(d['rank'] == 1).mean(), in_top_10=(d['rank'] <= 10).mean(),
                   median_rank=d['rank'].median(),
                   miss_among_10_nearest_all=(neighbour_all[cop[missed], winner[missed]] < 10).mean(),
                   miss_among_10_nearest_identified=(neighbour_identified[cop[missed], winner[missed]] < 10).mean())
        if kind in ('recovery', 'grid'):
            per_cop = d.groupby('cop').in_set.mean()
            row.update(in_set=d.in_set.mean(), mean_set_size=d.set_size.mean(),
                       cops_never_in_set=int((per_cop == 0).sum()), cops_in_set_majority=int((per_cop > .5).sum()))
        rows.append(row)

    table = pd.DataFrame(rows).sort_values(['rows', 'vector'])
    with open(f'{DIR}/report.txt', 'w') as f:
        f.write('One row per output vector and number of rows. smallest_loss: share of benchmarks whose own\n'
                'configuration has the smallest mean loss. in_set: share whose own configuration the\n'
                'confidence set keeps. miss_among_10_nearest_*: among benchmarks where another configuration\n'
                'has the smallest mean loss, the share where it is one of the 10 configurations nearest\n'
                'the right one, in all nine parameters scaled to the unit interval or in the seven\n'
                f'the responses identify ({", ".join(IDENTIFIED)}).\n\n')
        f.write(table.round(3).to_string(index=False))
        f.write('\n')


if __name__ == '__main__':
    main()
