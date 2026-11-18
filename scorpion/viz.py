"""Shared chart styling.

Colors come from a validated categorical palette: slots 1-3 clear the
colorblind-separation and normal-vision floors on this surface for every pair,
so no chart here uses more than three categorical series at once. Status colors
(pass / fail) are always paired with a label or icon, never carrying meaning on
hue alone.
"""

from __future__ import annotations

import plotly.graph_objects as go

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#86857f"
GRID = "#e8e7e2"

# categorical slots 1-3
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"

# status (fixed, never themed)
GOOD = "#0ca30c"
CRITICAL = "#d03b3b"

# one-hue sequential ramp, light -> dark, for ordered stages
BLUE_RAMP = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]

FONT = "ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif"


def style(fig: go.Figure, *, height: int = 360, showlegend: bool = False) -> go.Figure:
    """Apply the shared chrome: recessive grid and axes, quiet ink, no clutter."""
    fig.update_layout(
        height=height,
        showlegend=showlegend,
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, size=13, color=INK_SECONDARY),
        margin=dict(l=8, r=8, t=28, b=8),
        hoverlabel=dict(font=dict(family=FONT, size=12), bgcolor="#ffffff"),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            x=0,
            title_text="",
            font=dict(color=INK_SECONDARY),
        ),
        bargap=0.35,
    )
    fig.update_xaxes(
        showgrid=False,
        zeroline=False,
        linecolor=GRID,
        ticks="outside",
        tickcolor=GRID,
        title_font=dict(size=12, color=INK_MUTED),
    )
    fig.update_yaxes(
        gridcolor=GRID,
        zeroline=False,
        linecolor=GRID,
        ticks="",
        title_font=dict(size=12, color=INK_MUTED),
    )
    return fig


def value_labels(fig: go.Figure, template: str = "%{text}") -> go.Figure:
    """Direct labels on the data ends, in text ink rather than the series color."""
    fig.update_traces(
        textposition="outside",
        texttemplate=template,
        textfont=dict(color=INK, size=12, family=FONT),
        cliponaxis=False,
    )
    return fig


def rounded_bars(fig: go.Figure, radius: int = 4) -> go.Figure:
    """4px rounded data-ends, anchored to the baseline."""
    fig.update_traces(marker_cornerradius=radius, selector=dict(type="bar"))
    return fig
