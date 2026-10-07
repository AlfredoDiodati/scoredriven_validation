# How cop_0191 fares in the 1000 by 1000 Monte Carlo, against the whole grid,
# from the rows montecarlo/sweep_grid.c writes.
#
# Three tables per model, each giving the share of benchmarks whose confidence
# set holds the configuration, holds it alone, and gives it the smallest mean
# loss:
#   cop_0191 as the truth      its own 1000 benchmarks, the configuration is cop_0191
#   cop_0191 as a wrong answer the other 999,000 benchmarks, the configuration is
#                              cop_0191 although it did not generate the data
#   every true configuration   all 1,000,000 benchmarks, the configuration is the
#                              benchmark's own
#
# Reads montecarlo/out/sweep_grid.csv.gz, streamed row by row. Writes the three
# tables in percent to montecarlo/out/sweep_grid_cop_0191.txt and, as LaTeX
# tables, to montecarlo/out/sweep_grid_cop_0191.tex. Nothing printed.
#
# usage: python montecarlo/sweep_grid_cop_0191.py
import csv
import gzip
from collections import Counter, defaultdict

ROWS_PATH = 'montecarlo/out/sweep_grid.csv.gz'
REPORT_PATH = 'montecarlo/out/sweep_grid_cop_0191.txt'
LATEX_PATH = 'montecarlo/out/sweep_grid_cop_0191.tex'
COP = 191
MODEL_NAMES = {'qvarma': 'qvarma', 'lp_lin': 'linear', 'lp_s1': 'state 1', 'lp_s2': 'state 2',
               'lp_nl': 'both states'}
OUTCOMES = ('in_set', 'alone', 'smallest_loss')
OUTCOME_NAMES = ('in the set', 'alone in the set', 'smallest mean loss')


def tally(counts, members, lowest_loss_cop, cop):
    counts['benchmarks'] += 1
    counts['in_set'] += cop in members
    counts['alone'] += members == [cop]
    counts['smallest_loss'] += lowest_loss_cop == cop


def count_by_model():
    truth = defaultdict(Counter)
    wrong_answer = defaultdict(Counter)
    every_truth = defaultdict(Counter)
    with gzip.open(ROWS_PATH, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            model = row['model']
            benchmark_cop = int(row['benchmark_cop'])
            members = [int(member) for member in row['set_members'].split()]
            lowest_loss_cop = int(row['lowest_loss_cop'])
            tally(every_truth[model], members, lowest_loss_cop, benchmark_cop)
            if benchmark_cop == COP:
                tally(truth[model], members, lowest_loss_cop, COP)
            else:
                tally(wrong_answer[model], members, lowest_loss_cop, COP)
    return truth, wrong_answer, every_truth


def percent_rows(by_model):
    """One row per model: its name and the three shares, in percent, two decimals."""
    rows = []
    for model, counts in by_model.items():
        n = counts['benchmarks']
        rows.append((MODEL_NAMES[model], [f'{100 * counts[k] / n:.2f}' for k in OUTCOMES]))
    return rows


def write_text_table(report, title, by_model):
    report.write(f'{title}\n')
    report.write(f"order: {', '.join(OUTCOME_NAMES)}\n")
    for name, shares in percent_rows(by_model):
        report.write(f"{name}: {', '.join(share + '%' for share in shares)}\n")
    report.write('\n')


def write_latex_table(latex, caption, label, by_model):
    latex.write('\\begin{table}[htbp]\n\\centering\n')
    latex.write(f'\\caption{{{caption}}}\n\\label{{{label}}}\n')
    latex.write('\\begin{tabular}{lrrr}\n\\hline\n')
    latex.write('model & ' + ' & '.join(OUTCOME_NAMES) + ' \\\\\n\\hline\n')
    for name, shares in percent_rows(by_model):
        latex.write(name + ' & ' + ' & '.join(share + '\\%' for share in shares) + ' \\\\\n')
    latex.write('\\hline\n\\end{tabular}\n\\end{table}\n\n')


def main():
    truth, wrong_answer, every_truth = count_by_model()
    tables = [
        (truth, f'cop_{COP:04d} is the true configuration: share of its own '
         f'{truth["qvarma"]["benchmarks"]:,} benchmarks', 'truth'),
        (wrong_answer, f'cop_{COP:04d} is not the true configuration: share of the other '
         f'{wrong_answer["qvarma"]["benchmarks"]:,} benchmarks where cop_{COP:04d} is', 'wrong_answer'),
        (every_truth, f'every true configuration: share of all '
         f'{every_truth["qvarma"]["benchmarks"]:,} benchmarks where the true configuration is', 'every_truth'),
    ]
    with open(REPORT_PATH, 'w') as report:
        for by_model, title, _ in tables:
            write_text_table(report, title, by_model)
    with open(LATEX_PATH, 'w') as latex:
        for by_model, title, label in tables:
            caption = title.replace('_', '\\_').replace('%', '\\%')
            write_latex_table(latex, caption, f'tab:cop_{COP:04d}_{label}', by_model)


if __name__ == '__main__':
    main()
