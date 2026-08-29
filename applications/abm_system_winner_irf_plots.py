"""
Figures for out/abm_system_winner_irf.csv, which
applications/abm_system_winner_irf.c writes: the impulse responses of the
model the Model Confidence Set kept.

The same responses are drawn twice, under two identifications, because the
two answer different questions and neither one's shock labels mean what the
other's do.

out/abm_system_winner_irf_plots/sign_restricted/ is section 4.3 of Blazsek,
Escribano and Licht (2023). The structural shocks are whatever the rotations
satisfying the paper's Table 1 make them, so a shock is named for the sign
pattern that defines it and not for a variable: a rotation mixes the
reduced-form innovations, and column b of the impact matrix is the shock
whose impact signs match column b of Table 1. Two of the five shocks carry no
restriction and so carry no name either. Each column header states the
restriction on named variables that defines it, so nothing about the naming
is left implicit. The line is the 50th percentile over accepted rotations and
the band is the 10th and 90th.

out/abm_system_winner_irf_plots/recursive/ is the unrotated response, Q = I,
the one point in the identified set that Omega_inv's Cholesky orientation
picks out. Omega_inv is lower triangular, so shock b there is the part of
variable b's innovation orthogonal to the variables before it in the row
order the data already has - which is what makes every column name a
variable. It buys that at two prices, both stated on the figures: the answer
depends on the ordering, and there is no band, a single orientation having no
percentiles to take.

Per set:

  total_grid      every response to every shock, the whole matrix
  shock_<name>    one shock, its five responses at readable size
  cumulative      the response of the level, for the four series the data
                  carries as a difference
  decomposition   which channel carries the response - the transitory
                  I(0) recursion or the co-integrated random walk
  impact          the horizon-zero response matrix, where the sign
                  restrictions bind

and, in the sign-restricted set alone:

  identification  how much Table 1 actually pins down, the unrotated
                  response against the median over accepted rotations

Files are named for what they show, which the directory above them already
qualifies - so total_grid.pdf rather than
abm_system_winner_irf_plots_sign_restricted_total_grid.pdf. A panel is
addressed by the shock's column index rather than by the CSV's own
shock_name, since that column names the sign-restricted reading of the
column and means nothing under the recursive one. Nothing is printed.
"""

from collections import namedtuple
from pathlib import Path

import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

IRF_PATH = Path("out/abm_system_winner_irf.csv")
OUT_DIR = Path("out/abm_system_winner_irf_plots")

# Slots 1 to 3 of the data-visualization reference palette, the subset
# documented as clearing the colorblind and normal-vision separation floors
# over all pairs rather than only over adjacent ones.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
BAND_FILL = "rgba(42, 120, 214, 0.16)"

INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

SERIES_ORDER = [
    "GDP growth",
    "Energy demand growth",
    "Employment change",
    "Inflation",
    "Interest rate",
]

# What cumulating a response gives, for the four series the data carries as a
# first difference. The interest rate is a level already, so the sum of its
# responses is the response of nothing and it is left out of that figure.
CUMULATED_LEVEL = {
    "GDP growth": "GDP level",
    "Energy demand growth": "Energy demand level",
    "Employment change": "Employment level",
    "Inflation": "Price level",
}


def header(name, qualifier):
    """A column header on two lines: what the shock is called, and under it,
    smaller, the thing that makes the name mean something."""
    return (f"{name}<br><span style=\"font-size:10px;color:{INK_SOFT}\">"
            f"{qualifier}</span>")


# Table 1 of the paper, spelled out on the variables it restricts. Shocks 4
# and 5 exist because the system has five series; nothing restricts them, so
# nothing names them either.
SIGN_RESTRICTED_HEADERS = [
    header("Supply shock", "GDP +, inflation -"),
    header("Demand shock", "GDP +, inflation +, interest rate +"),
    header("Monetary policy shock", "GDP -, inflation -, interest rate +"),
    header("Unrestricted shock 4", "no restriction imposed"),
    header("Unrestricted shock 5", "no restriction imposed"),
]
SIGN_RESTRICTED_PHRASES = [
    "a supply shock",
    "a demand shock",
    "a monetary policy shock",
    "the fourth, unrestricted shock",
    "the fifth, unrestricted shock",
]

# Under the Cholesky orientation column b is variable b's own orthogonalized
# innovation, so the position in the ordering is what the name rests on and
# is stated beside it.
RECURSIVE_POSITION = ["first", "second", "third", "fourth", "fifth"]
# Written out rather than lower-cased from SERIES_ORDER, which would turn GDP
# into gDP and leave the interest rate without its article.
SERIES_IN_SENTENCE = [
    "GDP growth",
    "energy demand growth",
    "employment change",
    "inflation",
    "the interest rate",
]
RECURSIVE_HEADERS = [
    header(f"Shock to {name}", f"{position} in the recursive order")
    for name, position in zip(SERIES_IN_SENTENCE, RECURSIVE_POSITION)
]
RECURSIVE_PHRASES = [
    "a GDP growth shock",
    "an energy demand growth shock",
    "an employment change shock",
    "an inflation shock",
    "an interest rate shock",
]
RECURSIVE_FILE_STEMS = [
    "shock_gdp_growth",
    "shock_energy_demand_growth",
    "shock_employment_change",
    "shock_inflation",
    "shock_interest_rate",
]

Orientation = namedtuple(
    "Orientation",
    "directory headers phrases matrix_shocks panel_shocks file_stems banded line_label",
)

SIGN_RESTRICTED = Orientation(
    directory=OUT_DIR / "sign_restricted",
    headers=SIGN_RESTRICTED_HEADERS,
    phrases=SIGN_RESTRICTED_PHRASES,
    # The two unrestricted shocks are in the full matrix, where their width is
    # the point, and out of the figures that read one shock at a time.
    matrix_shocks=[1, 2, 3],
    panel_shocks=[1, 2, 3],
    file_stems=["shock_supply", "shock_demand", "shock_monetary"],
    banded=True,
    line_label="Median response over accepted rotations",
)

RECURSIVE = Orientation(
    directory=OUT_DIR / "recursive",
    headers=RECURSIVE_HEADERS,
    phrases=RECURSIVE_PHRASES,
    matrix_shocks=[1, 2, 3, 4, 5],
    panel_shocks=[1, 2, 3, 4, 5],
    file_stems=RECURSIVE_FILE_STEMS,
    banded=False,
    line_label="Response under the recursive (Cholesky) orientation",
)


def load_curves():
    """Every (component, shock index, response) as horizon-ordered lists, so
    a panel is one dictionary lookup rather than a filter over the frame."""
    frame = pl.read_csv(IRF_PATH)
    curves = {}
    for key, group in frame.group_by(["component", "shock", "response_name"]):
        ordered = group.sort("horizon")
        curves[tuple(key)] = {
            "horizon": ordered["horizon"].to_list(),
            "unrotated": ordered["unrotated"].to_list(),
            "lower": ordered["lower"].to_list(),
            "median": ordered["median"].to_list(),
            "upper": ordered["upper"].to_list(),
        }
    return curves


def line_of(curve, orientation):
    return curve["median"] if orientation.banded else curve["unrotated"]


def spread_of(curve, orientation):
    """What the y range has to cover: the band where there is one, the line
    itself where there is not."""
    if orientation.banded:
        return curve["lower"] + curve["upper"]
    return curve["unrotated"]


def style_axes(figure, x_title_rows=(), n_cols=1, n_rows=1):
    """Recessive grid and axes throughout, with the horizon named only on
    the bottom row so the label is not repeated once per panel."""
    figure.update_xaxes(
        showgrid=False,
        zeroline=False,
        linecolor=GRID,
        ticks="outside",
        tickcolor=GRID,
        ticklen=4,
        tickfont=dict(size=10, color=INK_SOFT),
        dtick=5,
    )
    figure.update_yaxes(
        gridcolor=GRID,
        zeroline=False,
        linecolor=GRID,
        ticks="outside",
        tickcolor=GRID,
        ticklen=4,
        tickfont=dict(size=10, color=INK_SOFT),
    )
    for row in x_title_rows:
        for col in range(1, n_cols + 1):
            figure.update_xaxes(
                title_text="Horizon (periods)",
                title_font=dict(size=11, color=INK_SOFT),
                row=row,
                col=col,
            )
    figure.update_layout(
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, size=11, color=INK),
        margin=dict(l=90, r=30, t=60, b=80),
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.05 if n_rows > 1 else -0.18,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            borderwidth=0,
            font=dict(size=11, color=INK_SOFT),
        ),
        hovermode=False,
    )


def add_zero_line(figure, horizon, row, col):
    figure.add_trace(
        go.Scatter(
            x=[horizon[0], horizon[-1]],
            y=[0, 0],
            mode="lines",
            line=dict(color=INK_SOFT, width=1, dash="dot"),
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def add_response(figure, curve, orientation, row, col, show_legend,
                 line_label=None):
    """The response, with its band underneath when the identification has
    one. The band's upper edge goes first and the lower edge fills back up to
    it, the order plotly's tonexty needs to shade the interval between."""
    if orientation.banded:
        figure.add_trace(
            go.Scatter(
                x=curve["horizon"],
                y=curve["upper"],
                mode="lines",
                line=dict(width=0, color=BAND_FILL),
                showlegend=False,
                legendgroup="band",
            ),
            row=row,
            col=col,
        )
        figure.add_trace(
            go.Scatter(
                x=curve["horizon"],
                y=curve["lower"],
                mode="lines",
                line=dict(width=0, color=BAND_FILL),
                fill="tonexty",
                fillcolor=BAND_FILL,
                name="10th to 90th percentile over accepted rotations",
                legendgroup="band",
                showlegend=show_legend,
            ),
            row=row,
            col=col,
        )
    figure.add_trace(
        go.Scatter(
            x=curve["horizon"],
            y=line_of(curve, orientation),
            mode="lines",
            line=dict(color=BLUE, width=2),
            name=line_label or orientation.line_label,
            legendgroup="line",
            showlegend=show_legend,
        ),
        row=row,
        col=col,
    )


def set_row_ranges(figure, row_values):
    """An explicit padded y range per row. Plotly's own autorange leaves a
    spike that reaches the top of the data sitting on the axis line, and a
    row here shares one axis across its columns, so the widest panel decides
    where that line is."""
    for row, values in enumerate(row_values, start=1):
        low, high = min(values), max(values)
        span = high - low
        pad = 0.08 * span if span > 0 else 1.0
        figure.update_yaxes(range=[low - pad, high + pad], row=row, col=1)


def add_matrix_headers(figure, columns, rows):
    """Column and row names once each, at the edges of the figure, instead
    of a repeated sentence in all n x m panel titles."""
    n_cols, n_rows = len(columns), len(rows)
    for index, name in enumerate(columns):
        figure.add_annotation(
            text=name,
            xref="paper",
            yref="paper",
            x=(index + 0.5) / n_cols,
            y=1.035,
            showarrow=False,
            font=dict(size=12, color=INK),
            xanchor="center",
            yanchor="bottom",
        )
    for index, name in enumerate(rows):
        figure.add_annotation(
            text=name,
            xref="paper",
            yref="paper",
            x=-0.055,
            y=1 - (index + 0.5) / n_rows,
            showarrow=False,
            textangle=-90,
            font=dict(size=12, color=INK),
            xanchor="center",
            yanchor="middle",
        )


def figure_total_grid(curves, orientation):
    """Every response to every shock, both identified and not. Rows share a
    y-axis, so the same variable's response is on one scale across the shocks
    and the panels can be read against each other."""
    shocks = list(range(1, len(orientation.headers) + 1))
    figure = make_subplots(
        rows=len(SERIES_ORDER),
        cols=len(shocks),
        shared_yaxes="rows",
        vertical_spacing=0.045,
        horizontal_spacing=0.025,
    )
    row_values = []
    for row, series in enumerate(SERIES_ORDER, start=1):
        spread = []
        for col, shock in enumerate(shocks, start=1):
            curve = curves[("total", shock, series)]
            add_zero_line(figure, curve["horizon"], row, col)
            add_response(figure, curve, orientation, row, col,
                         show_legend=(row == 1 and col == 1))
            spread += spread_of(curve, orientation)
        row_values.append(spread)
    style_axes(figure, x_title_rows=(len(SERIES_ORDER),),
               n_cols=len(shocks), n_rows=len(SERIES_ORDER))
    set_row_ranges(figure, row_values)
    add_matrix_headers(figure, orientation.headers, SERIES_ORDER)
    figure.update_layout(width=1250, height=1080, margin=dict(l=110, r=30, t=80, b=90))
    figure.write_image(orientation.directory / "total_grid.pdf")


def figure_one_shock(curves, orientation, shock, filename):
    """One shock at a size where the shape of each response is actually
    legible, five panels over two rows."""
    phrase = orientation.phrases[shock - 1]
    titles = [f"{series} to {phrase}" for series in SERIES_ORDER] + [""]
    figure = make_subplots(
        rows=2,
        cols=3,
        subplot_titles=titles,
        vertical_spacing=0.14,
        horizontal_spacing=0.075,
    )
    for index, series in enumerate(SERIES_ORDER):
        row, col = index // 3 + 1, index % 3 + 1
        curve = curves[("total", shock, series)]
        add_zero_line(figure, curve["horizon"], row, col)
        add_response(figure, curve, orientation, row, col, show_legend=(index == 0))
    figure.update_xaxes(visible=False, row=2, col=3)
    figure.update_yaxes(visible=False, row=2, col=3)
    style_axes(figure, n_cols=3, n_rows=2)
    # The third panel of the top row has no panel under it, so it carries its
    # own horizon label rather than borrowing the one below.
    for row, col in [(2, 1), (2, 2), (1, 3)]:
        figure.update_xaxes(title_text="Horizon (periods)",
                            title_font=dict(size=11, color=INK_SOFT), row=row, col=col)
    figure.update_yaxes(title_text="Percentage points",
                        title_font=dict(size=11, color=INK_SOFT), col=1)
    for annotation in figure.layout.annotations:
        annotation.font = dict(size=12, color=INK)
    figure.update_layout(width=1000, height=640,
                         margin=dict(l=80, r=30, t=50, b=120),
                         legend=dict(y=-0.14))
    figure.write_image(orientation.directory / filename)


def figure_cumulative(curves, orientation):
    """The response of the level rather than of the change. Four of the five
    series are first differences, so what a reader wants to know about them
    is where the level ends up, which is the response summed over the
    horizon. The band is the band of the cumulated path itself, taken over
    the accepted rotations by qvarma.h rather than summed from the reported
    band here - a quantile of a sum is not the sum of quantiles."""
    series_order = list(CUMULATED_LEVEL)
    shocks = orientation.matrix_shocks
    figure = make_subplots(
        rows=len(series_order),
        cols=len(shocks),
        shared_yaxes="rows",
        vertical_spacing=0.055,
        horizontal_spacing=0.04,
    )
    row_values = []
    for row, series in enumerate(series_order, start=1):
        spread = []
        for col, shock in enumerate(shocks, start=1):
            curve = curves[("cumulative", shock, series)]
            spread += spread_of(curve, orientation)
            add_zero_line(figure, curve["horizon"], row, col)
            add_response(figure, curve, orientation, row, col,
                         show_legend=(row == 1 and col == 1),
                         line_label="Median cumulated response"
                         if orientation.banded else None)
        row_values.append(spread)
    style_axes(figure, x_title_rows=(len(series_order),),
               n_cols=len(shocks), n_rows=len(series_order))
    set_row_ranges(figure, row_values)
    add_matrix_headers(figure, [orientation.headers[s - 1] for s in shocks],
                       [CUMULATED_LEVEL[series] for series in series_order])
    figure.update_layout(width=180 + 270 * len(shocks), height=900,
                         margin=dict(l=110, r=30, t=80, b=90))
    figure.write_image(orientation.directory / "cumulative.pdf")


def figure_decomposition(curves, orientation):
    """Where a response comes from. Under any one rotation the total is the
    impact response plus the transitory I(0) recursion plus the co-integrated
    random walk, and the last two are zero on different blocks of the system,
    so this says which of the two channels each variable rides on. In the
    sign-restricted set each line is that component's own median over the
    accepted rotations, so the three need not add up exactly - a median is
    not additive. The total is drawn widest and first, since a component that
    carries the whole response sits exactly on top of it."""
    shocks = orientation.matrix_shocks
    figure = make_subplots(
        rows=len(SERIES_ORDER),
        cols=len(shocks),
        shared_yaxes="rows",
        vertical_spacing=0.045,
        horizontal_spacing=0.04,
    )
    parts = [
        ("total", BLUE, "Total response", 4),
        ("stationary", ORANGE, "Transitory component", 1.8),
        ("cointegrated", AQUA, "Co-integrated component", 1.8),
    ]
    row_values = []
    for row, series in enumerate(SERIES_ORDER, start=1):
        spread = []
        for col, shock in enumerate(shocks, start=1):
            add_zero_line(figure, curves[("total", shock, series)]["horizon"], row, col)
            for component, color, label, width in parts:
                curve = curves[(component, shock, series)]
                values = line_of(curve, orientation)
                spread += values
                figure.add_trace(
                    go.Scatter(
                        x=curve["horizon"],
                        y=values,
                        mode="lines",
                        line=dict(color=color, width=width),
                        opacity=0.35 if component == "total" else 1.0,
                        name=label,
                        legendgroup=label,
                        showlegend=(row == 1 and col == 1),
                    ),
                    row=row,
                    col=col,
                )
        row_values.append(spread)
    style_axes(figure, x_title_rows=(len(SERIES_ORDER),),
               n_cols=len(shocks), n_rows=len(SERIES_ORDER))
    set_row_ranges(figure, row_values)
    add_matrix_headers(figure, [orientation.headers[s - 1] for s in shocks], SERIES_ORDER)
    figure.update_layout(width=180 + 270 * len(shocks), height=1080,
                         margin=dict(l=110, r=30, t=80, b=90))
    figure.write_image(orientation.directory / "decomposition.pdf")


def figure_identification(curves, orientation):
    """What the sign restrictions buy. The unrotated response is one point in
    the identified set and the band is the part of that set Table 1 keeps, so
    the distance between the two is how far the restrictions move the answer
    away from the Cholesky orientation - and it is also the distance between
    this set of figures and the recursive one beside it."""
    shocks = orientation.matrix_shocks
    figure = make_subplots(
        rows=len(SERIES_ORDER),
        cols=len(shocks),
        shared_yaxes="rows",
        vertical_spacing=0.045,
        horizontal_spacing=0.04,
    )
    row_values = []
    for row, series in enumerate(SERIES_ORDER, start=1):
        spread = []
        for col, shock in enumerate(shocks, start=1):
            curve = curves[("total", shock, series)]
            first = row == 1 and col == 1
            spread += curve["lower"] + curve["upper"] + curve["unrotated"]
            add_zero_line(figure, curve["horizon"], row, col)
            add_response(figure, curve, orientation, row, col, show_legend=first)
            figure.add_trace(
                go.Scatter(
                    x=curve["horizon"],
                    y=curve["unrotated"],
                    mode="lines",
                    line=dict(color=ORANGE, width=2, dash="dash"),
                    name="Unrotated response, Cholesky orientation",
                    legendgroup="unrotated",
                    showlegend=first,
                ),
                row=row,
                col=col,
            )
        row_values.append(spread)
    style_axes(figure, x_title_rows=(len(SERIES_ORDER),),
               n_cols=len(shocks), n_rows=len(SERIES_ORDER))
    set_row_ranges(figure, row_values)
    add_matrix_headers(figure, [orientation.headers[s - 1] for s in shocks], SERIES_ORDER)
    figure.update_layout(width=180 + 270 * len(shocks), height=1080,
                         margin=dict(l=110, r=30, t=80, b=90))
    figure.write_image(orientation.directory / "identification.pdf")


def figure_impact(curves, orientation):
    """The horizon-zero response matrix, which is the only thing Table 1
    restricts directly. One x-axis for all the panels, since every entry is
    in the same units and the comparison across shocks is the point."""
    shocks = list(range(1, len(orientation.headers) + 1))
    figure = make_subplots(
        rows=1,
        cols=len(shocks),
        shared_yaxes=True,
        subplot_titles=orientation.headers,
        horizontal_spacing=0.02,
    )
    positions = list(range(len(SERIES_ORDER)))[::-1]
    for col, shock in enumerate(shocks, start=1):
        for position, series in zip(positions, SERIES_ORDER):
            curve = curves[("contemporaneous", shock, series)]
            first = col == 1 and series == SERIES_ORDER[0]
            if orientation.banded:
                figure.add_trace(
                    go.Scatter(
                        x=[curve["lower"][0], curve["upper"][0]],
                        y=[position, position],
                        mode="lines",
                        line=dict(color=BLUE, width=6),
                        opacity=0.28,
                        name="10th to 90th percentile over accepted rotations",
                        legendgroup="band",
                        showlegend=first,
                    ),
                    row=1,
                    col=col,
                )
                figure.add_trace(
                    go.Scatter(
                        x=[curve["median"][0]],
                        y=[position],
                        mode="markers",
                        marker=dict(color=BLUE, size=9),
                        name="Median impact",
                        legendgroup="median",
                        showlegend=first,
                    ),
                    row=1,
                    col=col,
                )
            figure.add_trace(
                go.Scatter(
                    x=[curve["unrotated"][0]],
                    y=[position],
                    mode="markers",
                    marker=dict(color=ORANGE if orientation.banded else BLUE, size=9,
                                symbol="diamond-open" if orientation.banded else "circle",
                                line=dict(width=2, color=ORANGE)),
                    name="Unrotated impact, Cholesky orientation"
                    if orientation.banded else orientation.line_label,
                    legendgroup="unrotated",
                    showlegend=first,
                ),
                row=1,
                col=col,
            )
        figure.add_vline(x=0, line=dict(color=INK_SOFT, width=1, dash="dot"),
                         row=1, col=col)
    style_axes(figure, n_cols=len(shocks), n_rows=1)
    figure.update_yaxes(
        tickmode="array",
        tickvals=positions,
        ticktext=SERIES_ORDER,
        showgrid=False,
        tickfont=dict(size=11, color=INK),
        range=[-0.6, len(SERIES_ORDER) - 0.4],
    )
    # One x range across the panels: every entry is in the same units, and how
    # large one shock's impact is next to another's is exactly what this figure
    # is for. Per-panel autoranging would rescale each until they all looked
    # the same size.
    figure.update_xaxes(matches="x", showgrid=True, gridcolor=GRID, dtick=None,
                        title_text=None)
    figure.update_xaxes(
        title_text="Impact response (percentage points)",
        title_font=dict(size=11, color=INK_SOFT),
        row=1,
        col=len(shocks) // 2 + 1,
    )
    for annotation in figure.layout.annotations:
        annotation.font = dict(size=12, color=INK)
    figure.update_layout(width=1250, height=400,
                         margin=dict(l=150, r=30, t=60, b=110))
    figure.write_image(orientation.directory / "impact.pdf")


def main():
    curves = load_curves()
    for orientation in (SIGN_RESTRICTED, RECURSIVE):
        orientation.directory.mkdir(parents=True, exist_ok=True)
        figure_total_grid(curves, orientation)
        for shock, stem in zip(orientation.panel_shocks, orientation.file_stems):
            figure_one_shock(curves, orientation, shock, f"{stem}.pdf")
        figure_cumulative(curves, orientation)
        figure_decomposition(curves, orientation)
        figure_impact(curves, orientation)
    figure_identification(curves, SIGN_RESTRICTED)


if __name__ == "__main__":
    main()
