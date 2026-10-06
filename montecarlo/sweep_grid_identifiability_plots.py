"""
Figures for the rows montecarlo/sweep_grid.c writes: how easy each configuration
is to identify when one of its own runs stands in for the real data.

Every configuration (CoP) is the benchmark in 1000 confidence sets, one per run.
Two bar charts per auxiliary model, one bar per CoP on the x axis, both on a
0 to 100 axis so they can be read against each other:

  <model>_true_cop           the share of that CoP's 1000 confidence sets that
                             contain it. High means the procedure finds the CoP
                             when it is the right answer, low means it does not

  <model>_top_other_cop      for the same CoP, the share of its 1000 confidence
                             sets that contain the one other CoP appearing in
                             them most often. Hovering a bar names that CoP. A bar
                             as tall as the first chart's says another CoP is
                             as easy to confuse with the right answer as the
                             right answer is to find

A dashed line marks the mean over CoPs in each figure. The CoP numbers the bars
are keyed by are the numbers in the CoP names, so cop_0191 is at x = 191.

The numbers behind the figures are written once to
montecarlo/out/sweep_grid_identifiability.csv, one row per model and CoP: the two shares
and the other CoP's number. A model is drawn when its rows are on disk, as
montecarlo/out/sweep_grid_<model>.csv.gz or inside the combined
montecarlo/out/sweep_grid.csv.gz that bin/sweep_grid writes when every model is
done. Nothing is printed.
"""

from pathlib import Path

import polars as pl
import plotly.graph_objects as go

OUT_DIR = Path("montecarlo/out")
FIGURE_DIR = OUT_DIR / "sweep_grid_identifiability_plots"
TABLE_PATH = OUT_DIR / "sweep_grid_identifiability.csv"
COMBINED_PATH = OUT_DIR / "sweep_grid.csv.gz"

MODELS = {
    "qvarma": "t-QVARMA",
    "lp_lin": "Local projection, linear",
    "lp_s1": "Local projection, state 1",
    "lp_s2": "Local projection, state 2",
    "lp_nl": "Local projection, both states",
}

# The neutrals of applications/abm_system_winner_irf_plots.py, so the figures
# read as one set.
AQUA = "#1baf7a"
BLUE = "#1f5ad6"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"


def rows_of(model):
    """The benchmark rows of one model, or None when none are on disk."""
    path = OUT_DIR / f"sweep_grid_{model}.csv.gz"
    if not path.exists() and COMBINED_PATH.exists():
        path = COMBINED_PATH
    if not path.exists():
        return None
    rows = pl.read_csv(
        path,
        columns=["model", "benchmark_cop", "in_set", "set_members"],
        schema_overrides={"set_members": pl.String},
    )
    rows = rows.filter(pl.col("model") == model)
    return rows if rows.height else None


def shares_by_cop(rows):
    """Per CoP: the percentage of its confidence sets containing it, and the
    other CoP found in most of them with that CoP's percentage. A tie goes to
    the lower CoP number."""
    runs = rows.group_by("benchmark_cop").agg(
        pl.len().alias("runs"),
        (pl.col("in_set").mean() * 100).alias("true_cop_percent"),
    )
    members = (
        rows.select("benchmark_cop", pl.col("set_members").str.split(" ").alias("member"))
        .explode("member", empty_as_null=True)
        .with_columns(pl.col("member").cast(pl.Int64))
        .filter(pl.col("member") != pl.col("benchmark_cop"))
    )
    top_other = (
        members.group_by("benchmark_cop", "member")
        .agg(pl.len().alias("appearances"))
        .sort(["benchmark_cop", "appearances", "member"], descending=[False, True, False])
        .group_by("benchmark_cop", maintain_order=True)
        .first()
    )
    return (
        runs.join(top_other, on="benchmark_cop", how="left")
        .with_columns(
            pl.col("member").alias("top_other_cop"),
            (pl.col("appearances").fill_null(0) / pl.col("runs") * 100).alias("top_other_percent"),
        )
        .select("benchmark_cop", "runs", "true_cop_percent", "top_other_cop", "top_other_percent")
        .sort("benchmark_cop")
    )


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


def figure_true_cop(table, model):
    figure = go.Figure(go.Bar(
        x=table.get_column("benchmark_cop").to_list(),
        y=table.get_column("true_cop_percent").to_list(),
        marker=dict(color=AQUA, line=dict(width=0)),
        hovertemplate="cop_%{x:04d}: %{y:.1f}%<extra></extra>",
    ))
    mean = table.get_column("true_cop_percent").mean()
    style_bars(figure, f"{MODELS[model]}: share of a CoP's own confidence sets that contain it",
               "% of confidence sets containing the true CoP", mean)
    figure.write_image(FIGURE_DIR / f"{model}_true_cop.pdf")


def figure_top_other_cop(table, model):
    figure = go.Figure(go.Bar(
        x=table.get_column("benchmark_cop").to_list(),
        y=table.get_column("top_other_percent").to_list(),
        marker=dict(color=BLUE, line=dict(width=0)),
        customdata=table.get_column("top_other_cop").to_list(),
        hovertemplate="benchmark cop_%{x:04d}, most frequent other CoP cop_%{customdata:04d}: %{y:.1f}%<extra></extra>",
    ))
    mean = table.get_column("top_other_percent").mean()
    style_bars(figure, f"{MODELS[model]}: share of a CoP's confidence sets that contain the most frequent other CoP",
               "% of confidence sets containing the most frequent other CoP", mean)
    figure.write_image(FIGURE_DIR / f"{model}_top_other_cop.pdf")


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    tables = []
    for model in MODELS:
        rows = rows_of(model)
        if rows is None:
            continue
        table = shares_by_cop(rows)
        figure_true_cop(table, model)
        figure_top_other_cop(table, model)
        tables.append(table.with_columns(pl.lit(model).alias("model")).select(
            "model", "benchmark_cop", "runs", "true_cop_percent", "top_other_cop", "top_other_percent"))
    pl.concat(tables).write_csv(TABLE_PATH, float_precision=3)


if __name__ == "__main__":
    main()
