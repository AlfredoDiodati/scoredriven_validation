"""
Figures for out/montecarlo/seri_against_naive_by_configuration_identifiability.csv, which
bin/seri_against_naive_by_configuration writes: how easy each configuration is
to identify when one of its own runs stands in for the real data, under the
Seri grid of montecarlo/seri_true_configuration_recovery.c and under the
discarded naive-loss grid of old/naive_mae_loss/montecarlo/sweep_grid.c.

The figures are those of old/naive_mae_loss/montecarlo/sweep_grid_identifiability_plots.py,
drawn for both grids with the same style so they can be read against each
other. Every configuration (CoP) is the benchmark in 1000 confidence sets, one
per run. Two bar charts per grid and auxiliary model, one bar per CoP on the
x axis, both on a 0 to 100 axis:

  <grid>_<model>_true_cop        the share of that CoP's 1000 confidence sets
                                 that contain it

  <grid>_<model>_top_other_cop   for the same CoP, the share of its 1000
                                 confidence sets that contain the one other CoP
                                 appearing in them most often

A dashed line marks the mean over CoPs in each figure. The CoP numbers the bars
are keyed by are the numbers in the CoP names, so cop_0191 is at x = 191.

Writes out/montecarlo/seri_against_naive_by_configuration_plots/. Nothing is printed.

usage: python montecarlo/seri_against_naive_by_configuration_plots.py
"""

from pathlib import Path

import polars as pl
import plotly.graph_objects as go

TABLE_PATH = Path("out/montecarlo/seri_against_naive_by_configuration_identifiability.csv")
FIGURE_DIR = Path("out/montecarlo/seri_against_naive_by_configuration_plots")

GRIDS = {
    "Seri": "Seri et al. loss and confidence set",
    "naive": "naive loss, bootstrap confidence set",
}

MODELS = {
    "qvarma": "t-QVARMA",
    "lp_lin": "Local projection, linear",
    "lp_s1": "Local projection, state 1",
    "lp_s2": "Local projection, state 2",
    "lp_nl": "Local projection, both states",
}

# The colours and neutrals of the naive grid's figures, so the two sets match.
AQUA = "#1baf7a"
BLUE = "#1f5ad6"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"


def style_bars(figure, title, y_title, mean):
    figure.add_hline(
        y=mean,
        line=dict(color=INK_SOFT, width=1.5, dash="dash"),
    )
    figure.update_xaxes(
        title_text="CoP",
        title_font=dict(size=11, color=INK_SOFT),
        showgrid=False,
        zeroline=False,
        linecolor=GRID,
        tickfont=dict(size=10, color=INK_SOFT),
        range=[0, 1001],
    )
    figure.update_yaxes(
        title_text=y_title,
        title_font=dict(size=11, color=INK_SOFT),
        gridcolor=GRID,
        zeroline=False,
        linecolor=GRID,
        tickfont=dict(size=10, color=INK_SOFT),
        range=[0, 100],
    )
    figure.update_layout(
        width=1400,
        height=520,
        title=dict(text=f"{title}. Dashed line: mean over CoPs, {mean:.1f}%", font=dict(size=14, color=INK), x=0.01),
        font=dict(family=FONT, color=INK),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        showlegend=False,
        bargap=0,
        margin=dict(l=70, r=30, t=70, b=60),
    )


def figure_true_cop(table, grid, model):
    figure = go.Figure(go.Bar(
        x=table.get_column("benchmark_cop").to_list(),
        y=table.get_column("true_cop_percent").to_list(),
        marker=dict(color=AQUA, line=dict(width=0)),
        hovertemplate="cop_%{x:04d}: %{y:.1f}%<extra></extra>",
    ))
    mean = table.get_column("true_cop_percent").mean()
    style_bars(figure, f"{MODELS[model]}, {GRIDS[grid]}: share of a CoP's own confidence sets that contain it",
               "% of confidence sets containing the true CoP", mean)
    figure.write_image(FIGURE_DIR / f"{grid}_{model}_true_cop.pdf")


def figure_top_other_cop(table, grid, model):
    figure = go.Figure(go.Bar(
        x=table.get_column("benchmark_cop").to_list(),
        y=table.get_column("top_other_percent").to_list(),
        marker=dict(color=BLUE, line=dict(width=0)),
        customdata=table.get_column("top_other_cop").to_list(),
        hovertemplate="benchmark cop_%{x:04d}, most frequent other CoP cop_%{customdata:04d}: %{y:.1f}%<extra></extra>",
    ))
    mean = table.get_column("top_other_percent").mean()
    style_bars(figure, f"{MODELS[model]}, {GRIDS[grid]}: share of a CoP's confidence sets that contain the most "
               "frequent other CoP", "% of confidence sets containing the most frequent other CoP", mean)
    figure.write_image(FIGURE_DIR / f"{grid}_{model}_top_other_cop.pdf")


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rows = pl.read_csv(TABLE_PATH)
    for grid in GRIDS:
        for model in MODELS:
            table = rows.filter((pl.col("grid") == grid) & (pl.col("model") == model)).sort("benchmark_cop")
            if not table.height:
                continue
            figure_true_cop(table, grid, model)
            figure_top_other_cop(table, grid, model)


if __name__ == "__main__":
    main()
