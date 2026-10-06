# How often each configuration turns up in confidence sets that are not its own,
# from the rows montecarlo/sweep_grid.c writes.
#
# Every configuration (CoP) is the benchmark in 1000 confidence sets, one per run.
# Per model and CoP, two rates:
#   own rate         the share of its own 1000 sets that contain it
#   off-target rate  the share of the other CoPs' 999,000 sets that contain it
# A CoP with a high off-target rate is one the procedure keeps when the right
# answer is something else.
#
# Reads montecarlo/out/sweep_grid.csv.gz, streamed row by row, so memory stays
# small. Writes montecarlo/out/sweep_grid_presence.txt, one block per model in
# the order the models appear in the file. Nothing printed.
#
# usage: python montecarlo/sweep_grid_presence.py
import csv
import gzip
from collections import defaultdict

ROWS_PATH = 'montecarlo/out/sweep_grid.csv.gz'
REPORT_PATH = 'montecarlo/out/sweep_grid_presence.txt'
TOP = 10


def count_by_model():
    """Per model: benchmarks per CoP, own-set hits per CoP, and appearances of
    each CoP in sets whose benchmark is another CoP."""
    benchmarks = defaultdict(lambda: defaultdict(int))
    own_hits = defaultdict(lambda: defaultdict(int))
    off_hits = defaultdict(lambda: defaultdict(int))
    with gzip.open(ROWS_PATH, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            model = row['model']
            cop = int(row['benchmark_cop'])
            benchmarks[model][cop] += 1
            own_hits[model][cop] += int(row['in_set'])
            for member in row['set_members'].split():
                if int(member) != cop:
                    off_hits[model][int(member)] += 1
    return benchmarks, own_hits, off_hits


def write_block(report, model, benchmarks, own_hits, off_hits):
    cops = sorted(benchmarks)
    total = sum(benchmarks.values())
    own = {c: own_hits[c] / benchmarks[c] for c in cops}
    off = {c: off_hits[c] / (total - benchmarks[c]) for c in cops}
    mean_own = sum(own.values()) / len(cops)
    mean_off = sum(off.values()) / len(cops)
    top = sorted(cops, key=lambda c: (-off[c], c))[:TOP]
    off_slots = sum(off_hits.values())

    report.write(f'{model}: mean own-CoP rate {mean_own:.4f}; mean off-target presence per CoP {mean_off:.4f}\n')
    report.write('  CoPs whose off-target presence rate exceeds their own recovery rate: '
                 f'{sum(off[c] > own[c] for c in cops)} of {len(cops)}\n')
    report.write(f'  CoPs whose off-target presence exceeds the mean own rate {mean_own:.4f}: '
                 f'{sum(off[c] > mean_own for c in cops)}\n')
    report.write(f'  top {TOP} by off-target presence (cop, off rate, own rate):\n')
    for c in top:
        report.write(f'    cop_{c:04d} {off[c]:.4f} {own[c]:.4f}\n')
    report.write(f'  share of all off-target set slots taken by the top {TOP}: '
                 f'{sum(off_hits[c] for c in top) / off_slots:.3f}\n')


def main():
    benchmarks, own_hits, off_hits = count_by_model()
    with open(REPORT_PATH, 'w') as report:
        for i, model in enumerate(benchmarks):
            if i:
                report.write('\n')
            write_block(report, model, benchmarks[model], own_hits[model], off_hits[model])


if __name__ == '__main__':
    main()
