"""
Figures for the tables tests/abm_system_winner_diagnostics.c writes.

Every histogram is over the fits of the configuration the Model Confidence Set
kept, one fit per replicate. Each bar is coloured by the share of its fits whose
optimizer reported convergence: yellow when all of them converged, blue when
none did, and a linear blend of the two in between. Each figure also marks the US benchmark,
the median over all fits and the parameter set abm_system_winner_irf averages
over the fits.

  parameters/<name> every fitted parameter, nu included, one figure each, on a
                    linear axis. The bars cover the fits between the 1st and
                    99th percentile and the figure says how many lie outside.
                    nu runs up to about 10^12 and its 99th percentile is near
                    1.5 million, so over that range most of its fits would fall
                    in the first bar; its bars stop at the 75th percentile
                    instead. When the US value falls outside the range the axis
                    is not stretched to reach it, which would squeeze the bars
                    to slivers; a note at the edge of the axis gives the value
                    instead

  parameters_overview
                    the same histograms side by side on one page, so the spread
                    of every parameter can be compared at a glance

  parameters_unconstrained/<name>, parameters_unconstrained_overview
                    the same two figures for theta, the unconstrained value the
                    optimizer steps, before the qvarma link maps it to the
                    parameter. The link only changes Phi_star (tanh), the
                    diagonal of Omega_inv (exp) and nu (exp plus two); for the
                    rest theta is the parameter. Every parameter uses the 1st to
                    99th percentile here, nu included, since log(nu - 2) has no
                    tail long enough to need a narrower range. The averaged
                    parameter set is marked at the mean of theta, the value the
                    link then maps

The script stops without drawing if the table names a different configuration
from the one out/abm_system_mcs_joint.csv keeps now. Nothing is printed.
"""

import math
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

PARAMETERS_PATH = Path("out/abm_system_winner_diagnostics_parameters.csv")
WINNER_PARAMETERS_PATH = Path("out/abm_system_winner_diagnostics_winner_parameters.csv")
WINNER_THETA_PATH = Path("out/abm_system_winner_diagnostics_winner_theta.csv")
MCS_PATH = Path("out/abm_system_mcs_joint.csv")
OUT_DIR = Path("out/abm_system_winner_diagnostics_plots")

# The neutrals of applications/abm_system_winner_irf_plots.py, so the two sets
# of figures read as one.
AQUA = "#1baf7a"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e6e5e1"
SURFACE = "#fcfcfb"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

# Ends of the bar colour scale. The yellow is darkened enough to stand out
# against the off-white background.
CONVERGED_YELLOW = (240, 190, 20)
NOT_CONVERGED_BLUE = (31, 90, 214)

# Range of the fits the per-parameter histograms bin, as quantiles, so a few
# extreme fits do not flatten every other bar into one.
PARAMETER_LOWER_QUANTILE = 0.01
PARAMETER_UPPER_QUANTILE = 0.99
PARAMETER_BINS = 50

OVERVIEW_COLUMNS = 6


@dataclass(frozen=True)
class Scale:
    """One scale the parameters are drawn on: where the fits and the reference
    values are read from and where the figures go."""
    fits_path: Path
    us_column: str
    averaged_column: str
    parameter_dir: str
    overview_name: str
    constrained: bool
    # Quantile range of the bars for parameters whose tail is too long for the
    # default one, keyed by parameter.
    quantile_range: dict = field(default_factory=dict)


CONSTRAINED = Scale(
    fits_path=WINNER_PARAMETERS_PATH,
    us_column="us_benchmark",
    averaged_column="averaged_unconstrained",
    parameter_dir="parameters",
    overview_name="parameters_overview.pdf",
    constrained=True,
    quantile_range={"nu": (0.01, 0.75)},
)

UNCONSTRAINED = Scale(
    fits_path=WINNER_THETA_PATH,
    us_column="us_benchmark_theta",
    averaged_column="averaged_theta",
    parameter_dir="parameters_unconstrained",
    overview_name="parameters_unconstrained_overview.pdf",
    constrained=False,
)


def check_winning_configuration(drawn):
    """The configuration the table belongs to, checked against the one the
    confidence set keeps now, since the two files are written by different
    runs."""
    kept = (
        pl.read_csv(MCS_PATH, columns=["model", "mean_loss", "in_set"], infer_schema_length=None)
        .filter(pl.col("in_set") == 1)
        .sort("mean_loss")
        .get_column("model")[0]
    )
    if not kept.startswith(f"{drawn}_"):
        raise SystemExit(
            f"{WINNER_PARAMETERS_PATH.parent} holds the fits of {drawn}, but {MCS_PATH} keeps {kept}: "
            f"run make study-abm_system_winner_diagnostics first"
        )


def load_winner(scale):
    """The winner's fits on one scale, one row per replicate with every
    parameter and the convergence flag, the reference table with the US and
    averaged values, the parameter names, and the configuration name."""
    reference = pl.read_csv(PARAMETERS_PATH)
    if not scale.fits_path.exists() or scale.averaged_column not in reference.columns:
        raise SystemExit(f"{scale.fits_path} or its columns in {PARAMETERS_PATH} are missing: "
                         f"run make study-abm_system_winner_diagnostics first")
    fits = pl.read_csv(scale.fits_path)
    configurations = fits.get_column("configuration").unique().to_list()
    if len(configurations) != 1:
        raise SystemExit(f"{scale.fits_path} holds more than one configuration")
    configuration = configurations[0]
    check_winning_configuration(configuration)

    parameters = [name for name in fits.columns if name not in ("configuration", "replicate", "is_converged")]
    missing = set(parameters) - set(reference.get_column("parameter").to_list())
    if missing:
        raise SystemExit(f"{PARAMETERS_PATH} has no row for {sorted(missing)}: the two tables come from different runs")
    return fits, reference, parameters, configuration


def reference_values(reference, parameter, scale):
    row = reference.filter(pl.col("parameter") == parameter)
    return row.get_column(scale.us_column)[0], row.get_column(scale.averaged_column)[0]


def link_of(reference, parameter):
    return reference.filter(pl.col("parameter") == parameter).get_column("link")[0]


def axis_title(reference, parameter, scale):
    if scale.constrained:
        return parameter
    return f"theta of {parameter}, unconstrained ({parameter} = {link_of(reference, parameter)})"


def overview_title(reference, parameter, scale):
    """The panel title, which names the link on the unconstrained scale when it
    is not the identity, and the quantile range when it is not the default one
    the legend states."""
    title = parameter
    if not scale.constrained and link_of(reference, parameter) != "theta":
        title = f"{parameter} = {link_of(reference, parameter)}"
    if parameter in scale.quantile_range:
        lower_quantile, upper_quantile = scale.quantile_range[parameter]
        title = f"{title}, {lower_quantile:.0%} to {upper_quantile:.0%} quantile"
    return title


def convergence_color(share):
    """Yellow for a bar whose fits all converged, blue for none, a linear blend
    between."""
    channels = [round(blue + share * (yellow - blue)) for yellow, blue in zip(CONVERGED_YELLOW, NOT_CONVERGED_BLUE)]
    return f"rgb({channels[0]},{channels[1]},{channels[2]})"


def convergence_colorscale():
    return [[0, convergence_color(0)], [0.5, convergence_color(0.5)], [1, convergence_color(1)]]


def binned_bars(values, converged, start, bin_size, n_bins):
    """Counts and converged share per bin for the values in [start, start +
    n_bins * bin_size], a value on the upper edge going to the last bin."""
    bins = (
        pl.DataFrame({"value": values, "converged": converged})
        .with_columns(((pl.col("value") - start) / bin_size).floor().cast(pl.Int64).clip(0, n_bins - 1).alias("bin"))
        .group_by("bin")
        .agg(pl.len().alias("count"), pl.col("converged").mean().alias("share"))
        .sort("bin")
    )
    centers = [start + (b + 0.5) * bin_size for b in bins.get_column("bin").to_list()]
    shares = bins.get_column("share").to_list()
    counts = bins.get_column("count").to_list()
    return centers, counts, shares


def bar_trace(centers, counts, shares, bin_size, name, showlegend=True):
    return go.Bar(
        x=centers,
        y=counts,
        width=bin_size,
        marker=dict(color=[convergence_color(share) for share in shares], line=dict(width=0)),
        name=name,
        showlegend=showlegend,
    )


def colorbar_trace(colorbar_y=0.5, length=0.8):
    """An invisible point carrying the colour scale, so the figure explains the
    bar colours."""
    return go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            color=[0],
            cmin=0,
            cmax=1,
            colorscale=convergence_colorscale(),
            showscale=True,
            colorbar=dict(
                title=dict(text="Share of converged<br>fits in the bar", font=dict(size=10, color=INK_SOFT)),
                tickvals=[0, 0.5, 1],
                ticktext=["0", "0.5", "1"],
                tickfont=dict(size=10, color=INK_SOFT),
                thickness=12,
                len=length,
                y=colorbar_y,
                outlinewidth=0,
            ),
        ),
        showlegend=False,
    )


def add_us_marker(figure, summary, label, font_size, **position):
    """The US benchmark as a vertical line when it lies on the axis, and as a
    note at the nearer edge of the axis when it does not."""
    low, high = summary["axis_range"]
    us_value = summary["us_value"]
    if low <= us_value <= high:
        add_marker_line(figure, us_value, INK, "solid", label, summary["top"], **position)
        return
    add_marker_line(figure, None, INK, "solid", label, summary["top"], **position)
    above = us_value > high
    # Subplot axes are numbered row first from 1, and the first is plain x and y.
    axis_number = (position["row"] - 1) * OVERVIEW_COLUMNS + position["col"] if position else 1
    axis_suffix = "" if axis_number == 1 else str(axis_number)
    figure.add_annotation(
        x=high if above else low,
        y=summary["top"] * 0.95,
        xref=f"x{axis_suffix}",
        yref=f"y{axis_suffix}",
        text=f"US {us_value:.3g} &#8594;" if above else f"&#8592; US {us_value:.3g}",
        showarrow=False,
        xanchor="right" if above else "left",
        font=dict(size=font_size, color=INK),
    )


def add_marker_line(figure, x, color, dash, label, top, **position):
    """A vertical reference line drawn as a trace, so it gets a legend entry."""
    figure.add_trace(
        go.Scatter(
            x=[x, x],
            y=[0, top],
            mode="lines",
            line=dict(color=color, width=2, dash=dash),
            name=label,
            showlegend=bool(label),
        ),
        **position,
    )


def style_histogram(figure, x_title, top):
    figure.update_xaxes(
        title_text=x_title,
        title_font=dict(size=11, color=INK_SOFT),
        showgrid=False,
        zeroline=False,
        linecolor=GRID,
        ticks="outside",
        tickcolor=GRID,
        ticklen=4,
        tickfont=dict(size=10, color=INK_SOFT),
    )
    figure.update_yaxes(
        title_text="Number of fits",
        title_font=dict(size=11, color=INK_SOFT),
        range=[0, top],
        gridcolor=GRID,
        zeroline=False,
        linecolor=GRID,
        ticks="outside",
        tickcolor=GRID,
        ticklen=4,
        tickfont=dict(size=10, color=INK_SOFT),
    )
    figure.update_layout(
        width=1000,
        height=520,
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, size=11, color=INK),
        margin=dict(l=80, r=120, t=40, b=150),
        bargap=0,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            borderwidth=0,
            font=dict(size=11, color=INK_SOFT),
        ),
        hovermode=False,
    )


def converged_count_label(fits):
    return f"{fits.get_column('is_converged').sum()} of {fits.height} converged"


def file_stem(parameter):
    """Psi_star1[0,2] -> Psi_star1_0_2."""
    return parameter.replace("[", "_").replace("]", "").replace(",", "_")


def parameter_summary(fits, reference, parameter, scale):
    """The bars to draw for one parameter, the axis range and the three
    reference values."""
    values = fits.get_column(parameter)
    us_value, averaged_value = reference_values(reference, parameter, scale)
    median_value = values.median()

    lower_quantile, upper_quantile = scale.quantile_range.get(
        parameter, (PARAMETER_LOWER_QUANTILE, PARAMETER_UPPER_QUANTILE))
    low = values.quantile(lower_quantile, interpolation="linear")
    high = values.quantile(upper_quantile, interpolation="linear")
    if high <= low:
        spread = abs(low) * 1e-3 if low != 0 else 1e-3
        low, high = low - spread, high + spread
    shown = fits.filter((pl.col(parameter) >= low) & (pl.col(parameter) <= high))
    bin_size = (high - low) / PARAMETER_BINS
    centers, counts, shares = binned_bars(shown.get_column(parameter), shown.get_column("is_converged"),
                                          low, bin_size, PARAMETER_BINS)

    return dict(
        centers=centers,
        counts=counts,
        shares=shares,
        lower_quantile=lower_quantile,
        upper_quantile=upper_quantile,
        n_shown=shown.height,
        n_outside=fits.height - shown.height,
        bin_size=bin_size,
        top=max(counts) * 1.08,
        axis_range=[min(low, averaged_value) - bin_size, max(high, averaged_value) + bin_size],
        us_value=us_value,
        averaged_value=averaged_value,
        median_value=median_value,
    )


def figure_parameter_distributions(fits, reference, parameters, configuration, scale):
    parameter_dir = OUT_DIR / scale.parameter_dir
    parameter_dir.mkdir(parents=True, exist_ok=True)

    for parameter in parameters:
        summary = parameter_summary(fits, reference, parameter, scale)
        value_name = parameter if scale.constrained else f"theta of {parameter}"
        figure = go.Figure()
        figure.add_trace(bar_trace(
            summary["centers"], summary["counts"], summary["shares"], summary["bin_size"],
            f"Fits of {configuration} between the {summary['lower_quantile']:.0%} and "
            f"{summary['upper_quantile']:.0%} quantiles ({summary['n_shown']} of {fits.height}; "
            f"{summary['n_outside']} outside not shown; {converged_count_label(fits)})",
        ))
        figure.add_trace(colorbar_trace())
        add_us_marker(figure, summary, f"US benchmark, {value_name} = {summary['us_value']:.4g}", 11)
        add_marker_line(figure, summary["median_value"], INK_SOFT, "dot",
                        f"Median over all fits, {value_name} = {summary['median_value']:.4g}", summary["top"])
        add_marker_line(figure, summary["averaged_value"], AQUA, "dash",
                        f"Averaged parameter set used for the IRF, {value_name} = {summary['averaged_value']:.4g}",
                        summary["top"])

        style_histogram(figure, axis_title(reference, parameter, scale), summary["top"])
        figure.update_xaxes(range=summary["axis_range"])
        figure.write_image(parameter_dir / f"{file_stem(parameter)}.pdf")


def figure_parameters_overview(fits, reference, parameters, configuration, scale):
    rows = math.ceil(len(parameters) / OVERVIEW_COLUMNS)
    figure = make_subplots(
        rows=rows,
        cols=OVERVIEW_COLUMNS,
        subplot_titles=[overview_title(reference, parameter, scale) for parameter in parameters],
        horizontal_spacing=0.04,
        vertical_spacing=0.045,
    )

    for index, parameter in enumerate(parameters):
        position = dict(row=index // OVERVIEW_COLUMNS + 1, col=index % OVERVIEW_COLUMNS + 1)
        summary = parameter_summary(fits, reference, parameter, scale)
        first = index == 0
        figure.add_trace(
            bar_trace(
                summary["centers"], summary["counts"], summary["shares"], summary["bin_size"],
                f"Fits of {configuration}{'' if scale.constrained else ', unconstrained theta,'} between the "
                f"{PARAMETER_LOWER_QUANTILE:.0%} and {PARAMETER_UPPER_QUANTILE:.0%} quantiles, one per replicate "
                f"({converged_count_label(fits)})",
                showlegend=first,
            ),
            **position,
        )
        add_us_marker(figure, summary, "US benchmark" if first else "", 8, **position)
        add_marker_line(figure, summary["median_value"], INK_SOFT, "dot", "Median over all fits" if first else "",
                        summary["top"], **position)
        add_marker_line(figure, summary["averaged_value"], AQUA, "dash",
                        "Averaged parameter set used for the IRF" if first else "", summary["top"], **position)
        figure.update_xaxes(range=summary["axis_range"], **position)
        figure.update_yaxes(range=[0, summary["top"]], **position)

    figure.add_trace(colorbar_trace(colorbar_y=0.5, length=0.3), row=1, col=1)
    figure.update_traces(selector=dict(type="scatter", mode="lines"), line_width=1.5)
    figure.update_xaxes(showgrid=False, zeroline=False, linecolor=GRID, tickfont=dict(size=7, color=INK_SOFT), nticks=4)
    figure.update_yaxes(gridcolor=GRID, zeroline=False, linecolor=GRID, tickfont=dict(size=7, color=INK_SOFT), nticks=3)
    figure.update_annotations(font=dict(size=9, color=INK))
    figure.update_layout(
        width=1400,
        height=230 * rows + 120,
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, size=9, color=INK),
        margin=dict(l=40, r=120, t=40, b=90),
        bargap=0,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.02,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            borderwidth=0,
            font=dict(size=11, color=INK_SOFT),
        ),
        hovermode=False,
    )
    figure.write_image(OUT_DIR / scale.overview_name)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scales = (CONSTRAINED, UNCONSTRAINED)
    loaded = [load_winner(scale) for scale in scales]
    for scale, (fits, reference, parameters, configuration) in zip(scales, loaded):
        figure_parameters_overview(fits, reference, parameters, configuration, scale)
    for scale, (fits, reference, parameters, configuration) in zip(scales, loaded):
        figure_parameter_distributions(fits, reference, parameters, configuration, scale)


if __name__ == "__main__":
    main()
