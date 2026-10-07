"""
Histograms of the oracle fit metric v (docs/MONTECARLO_FIT_METRIC.md) from the
rows bin/fit_metric_oracle writes, two per model:

  <model>_all        every benchmark of the 1000 by 1000 Monte Carlo, every run
                     of every configuration, 1,000,000 values
  <model>_cop_0191   the 1000 benchmarks whose true configuration is cop_0191,
                     its runs 0 to 999

v lies between 0 and 1, so every histogram has the same 50 bins of width 0.02
on the same 0 to 1 axis, and the bar height is the percentage of that
histogram's benchmarks falling in the bin. A dashed line marks the median. The
counts are binned here rather than handed to plotly, so a figure carries 50
bars rather than a million points.

v = 1 / (1 + d / sigma) does not depend on the units of the responses, so the
histograms of different models can be read against each other.

Reads montecarlo/out/fit_metric_oracle.csv.gz. Writes one pdf per histogram to
montecarlo/out/fit_metric_v_histograms/. Needs polars, plotly and kaleido.
Nothing is printed.
"""

from pathlib import Path

import polars as pl
import plotly.graph_objects as go

OUT_DIR = Path("montecarlo/out")
ORACLE_PATH = OUT_DIR / "fit_metric_oracle.csv.gz"
FIGURE_DIR = OUT_DIR / "fit_metric_v_histograms"
COP = 191
BIN_WIDTH = 0.02
N_BINS = 50

MODELS = {
    "qvarma": "t-QVARMA",
    "lp_lin": "Local projection, linear",
    "lp_s1": "Local projection, state 1",
    "lp_s2": "Local projection, state 2",
    "lp_nl": "Local projection, both states",
}

# The palette of montecarlo/sweep_grid_identifiability_plots.py, so the figures
# read as one set.
AQUA = "#1baf7a"
BLUE = "#1f5ad6"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"


def percent_by_bin(v):
    """The percentage of the values in each of the N_BINS bins; v = 1 goes to
    the last bin."""
    counts = (
        pl.DataFrame({"v": v})
        .with_columns((pl.col("v") / BIN_WIDTH).floor().cast(pl.Int64).clip(0, N_BINS - 1).alias("bin"))
        .group_by("bin")
        .agg(pl.len().alias("count"))
    )
    percent = [0.0] * N_BINS
    for row in counts.iter_rows(named=True):
        percent[row["bin"]] = 100 * row["count"] / len(v)
    return percent


def histogram(v, title, color, path):
    centers = [(i + 0.5) * BIN_WIDTH for i in range(N_BINS)]
    median = v.median()
    figure = go.Figure(go.Bar(
        x=centers,
        y=percent_by_bin(v),
        width=BIN_WIDTH,
        marker=dict(color=color, line=dict(width=0)),
        hovertemplate="v %{x:.2f}: %{y:.2f}%<extra></extra>",
    ))
    figure.add_vline(x=median, line=dict(color=INK_SOFT, width=1.5, dash="dash"))
    figure.update_xaxes(
        title_text="v",
        title_font=dict(size=11, color=INK_SOFT),
        showgrid=False,
        zeroline=False,
        linecolor=GRID,
        tickfont=dict(size=10, color=INK_SOFT),
        range=[0, 1],
    )
    figure.update_yaxes(
        title_text="% of benchmarks",
        title_font=dict(size=11, color=INK_SOFT),
        gridcolor=GRID,
        zeroline=False,
        linecolor=GRID,
        tickfont=dict(size=10, color=INK_SOFT),
    )
    figure.update_layout(
        width=900,
        height=520,
        title=dict(text=f"{title}, {len(v):,} benchmarks. Dashed line: median, {median:.3f}",
                   font=dict(size=14, color=INK), x=0.01),
        font=dict(family=FONT, color=INK),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        showlegend=False,
        bargap=0,
        margin=dict(l=70, r=30, t=70, b=60),
    )
    figure.write_image(path)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rows = pl.read_csv(ORACLE_PATH, columns=["model", "benchmark_cop", "v"])
    for model, name in MODELS.items():
        of_model = rows.filter(pl.col("model") == model)
        histogram(of_model.get_column("v"), f"{name}: v over every benchmark", AQUA,
                  FIGURE_DIR / f"{model}_all.pdf")
        histogram(of_model.filter(pl.col("benchmark_cop") == COP).get_column("v"),
                  f"{name}: v over the runs of cop_{COP:04d}", BLUE, FIGURE_DIR / f"{model}_cop_{COP:04d}.pdf")


if __name__ == "__main__":
    main()
