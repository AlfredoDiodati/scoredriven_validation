# The raw-response baselines of docs/MONTECARLO_COMPRESSED_RESPONSE.md at full
# scale: the rows of montecarlo/out/sweep_grid.csv.gz whose benchmark is run 0 or
# run 1 of a configuration, in the columns montecarlo/compressed_response_recovery.c
# writes, so montecarlo/compressed_response_report.py reads them like any other
# recovery file. Nothing is recomputed; sweep_grid.c's own confidence sets are the
# baseline.
#
# Reads montecarlo/out/sweep_grid.csv.gz, streamed. Writes
# montecarlo/out/compressed_response/grid_<model>_raw_rows1000.csv, one per model
# in the file, rows in the order the file holds them. Nothing printed.
#
# usage: python montecarlo/compressed_response_grid_baseline.py
import csv
import gzip

ROWS_PATH = 'montecarlo/out/sweep_grid.csv.gz'
OUT_DIR = 'montecarlo/out/compressed_response'
BENCHMARK_RUNS = {'0', '1'}
COLUMNS = ['cop', 'run', 'in_set', 'rank', 'set_size', 'lowest_loss_cop', 'final_pvalue']


def main():
    outputs = {}
    with gzip.open(ROWS_PATH, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            if row['benchmark_run'] not in BENCHMARK_RUNS:
                continue
            model = row['model']
            if model not in outputs:
                handle = open(f'{OUT_DIR}/grid_{model}_raw_rows1000.csv', 'w', newline='')
                writer = csv.writer(handle, lineterminator='\n')
                writer.writerow(COLUMNS)
                outputs[model] = (handle, writer)
            # The p-value as a plain float, 0.388 rather than the file's 0.388000.
            outputs[model][1].writerow([row['benchmark_cop'], row['benchmark_run'], row['in_set'], row['rank'],
                                        row['set_size'], row['lowest_loss_cop'], repr(float(row['final_pvalue']))])
    for handle, _ in outputs.values():
        handle.close()


if __name__ == '__main__':
    main()
