"""
Figure for out/studies/response_variance_run_count_by_configuration.csv, which
bin/response_variance_run_count writes: how the summed response variance of a
configuration settles as runs are added from the last one backwards.

For configuration i and m runs (runs 1000-m .. 999), V_i(m) is the sum over
every entry of the response vector of the entry's sample variance, divisor
m - 1. The figure draws V_i(m) / V_i(1000) against m, one panel per model, as
bands across the 1000 configurations: smallest to largest, 5th to 95th
percentile, 25th to 75th percentile, and the median. Dashed lines at 0.95 and
1.05 mark the 5% tolerance of the study, a dotted line the 100 runs the
held-out-mean grid gives to the mean. Log scale on the ratio.

Writes out/studies/response_variance_run_count_plots/expanding_variance_ratio.pdf.
Nothing is printed.

usage: python studies/response_variance_run_count_plots.py
"""

from pathlib import Path

import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

TABLE_PATH = Path("out/studies/response_variance_run_count_by_configuration.csv")
FIGURE_DIR = Path("out/studies/response_variance_run_count_plots")

MODELS = {
    "lp_lin": "LP, linear",
    "lp_nl_causal_hp_estimated_lambda": "NLP, both states, causal HP cycle at the estimated lambda",
    "qvarma_growth_states": "t-QVARMA, both growth states",
}
MEAN_RUNS = 100
TOLERANCE = 0.05

BLUE = "#1f5ad6"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

# Lower and upper quantile of each band, widest first, and its opacity.
BANDS = [(0.0, 1.0, 0.12), (0.05, 0.95, 0.22), (0.25, 0.75, 0.38)]


def ratio_quantiles(table, model):
    """Quantiles across configurations of V_i(m) / V_i(n), one column per m."""
    rows = table.filter(pl.col("model") == model).drop("model", "cop")
    columns = rows.columns
    full = columns[-1]
    ratio = rows.select([pl.col(c) / pl.col(full) for c in columns])
    runs = [int(c.split("_")[1]) for c in columns]
    levels = sorted({q for band in BANDS for q in band[:2]} | {0.5})
    quantile = {q: ratio.select(pl.all().quantile(q, interpolation="linear")).row(0) for q in levels}
    return runs, quantile


def band(figure, row, runs, lower, upper, opacity, name):
    figure.add_trace(go.Scatter(x=runs, y=lower, mode="lines", line=dict(width=0), showlegend=False,
                                hoverinfo="skip"), row=row, col=1)
    figure.add_trace(go.Scatter(x=runs, y=upper, mode="lines", line=dict(width=0), fill="tonexty",
                                fillcolor=f"rgba(31, 90, 214, {opacity})", name=name, showlegend=row == 1,
                                legendgroup=name), row=row, col=1)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    table = pl.read_csv(TABLE_PATH)
    figure = make_subplots(rows=len(MODELS), cols=1, shared_xaxes=True, vertical_spacing=0.06,
                           subplot_titles=list(MODELS.values()))
    for row, model in enumerate(MODELS, start=1):
        runs, quantile = ratio_quantiles(table, model)
        for lower, upper, opacity in BANDS:
            name = "smallest to largest" if lower == 0 else f"{lower * 100:.0f}th to {upper * 100:.0f}th percentile"
            band(figure, row, runs, quantile[lower], quantile[upper], opacity, name)
        figure.add_trace(go.Scatter(x=runs, y=quantile[0.5], mode="lines", line=dict(color=BLUE, width=1.6),
                                    name="median", showlegend=row == 1, legendgroup="median"), row=row, col=1)
        for level in (1 - TOLERANCE, 1 + TOLERANCE):
            figure.add_hline(y=level, line=dict(color=INK_SOFT, width=1, dash="dash"), row=row, col=1)
        figure.add_vline(x=MEAN_RUNS, line=dict(color=INK_SOFT, width=1, dash="dot"), row=row, col=1)
        figure.update_yaxes(type="log", title_text="V(m) / V(1000)", title_font=dict(size=11, color=INK_SOFT),
                            gridcolor=GRID, zeroline=False, linecolor=GRID, tickfont=dict(size=10, color=INK_SOFT),
                            row=row, col=1)
    figure.update_xaxes(showgrid=False, zeroline=False, linecolor=GRID, tickfont=dict(size=10, color=INK_SOFT),
                        range=[2, 1000])
    figure.update_xaxes(title_text="m, runs counted from the last", title_font=dict(size=11, color=INK_SOFT),
                        row=len(MODELS), col=1)
    figure.update_annotations(font=dict(size=12, color=INK))
    figure.update_layout(
        width=1100,
        height=1150,
        title=dict(text="Summed response variance over the last m runs, relative to all 1000 runs, across the "
                        "1000 configurations. Dashed: 5% tolerance. Dotted: m = 100.",
                   font=dict(size=13, color=INK), x=0.01),
        font=dict(family=FONT, color=INK),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        legend=dict(orientation="h", y=-0.05, x=0, font=dict(size=11, color=INK_SOFT)),
        margin=dict(l=80, r=30, t=90, b=90),
    )
    figure.write_image(FIGURE_DIR / "expanding_variance_ratio.pdf")


if __name__ == "__main__":
    main()
