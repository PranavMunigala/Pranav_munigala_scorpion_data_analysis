"""Watch 665 attempted pairings narrow to 17 binders — and check the published labels."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from scorpion import data, viz

st.set_page_config(page_title="Binder criteria", page_icon="✅", layout="wide")
st.title("✅ From 665 attempts to 17 binders")
st.caption(
    "Every criterion re-derived from the raw per-replicate data, then checked "
    "against the labels published with the release."
)

p = data.pairs()
passes, published, disagreements = data.reproduction_check()

# --------------------------------------------------------------------------- #
# reproduction check
# --------------------------------------------------------------------------- #

m = st.columns(3)
m[0].metric("Pairs passing all four criteria", passes)
m[1].metric("Binders labelled in the release", published)
m[2].metric("Disagreements", disagreements)

if disagreements == 0:
    st.success(
        "The four criteria, applied to the raw data, reproduce the published set "
        "of binders exactly — same designs, same toxins, nothing extra and "
        "nothing missing.",
        icon="✅",
    )
else:
    st.error(f"{disagreements} pair(s) disagree with the published labels.", icon="⚠️")

st.divider()

# --------------------------------------------------------------------------- #
# the funnel
# --------------------------------------------------------------------------- #

st.subheader("Where the attempts fall away")

stages = [("All intended design–toxin pairs", np.ones(len(p), dtype=bool))]
mask = np.ones(len(p), dtype=bool)
for col, label in data.CRITERIA:
    mask = mask & p[col].to_numpy()
    stages.append((f"... and {label[0].lower() + label[1:]}", mask.copy()))

funnel = pd.DataFrame(
    {"stage": [s for s, _ in stages], "pairs": [int(m_.sum()) for _, m_ in stages]}
)

fig = go.Figure(
    go.Bar(
        x=funnel.pairs,
        y=funnel.stage,
        orientation="h",
        marker_color=viz.BLUE_RAMP,
        text=funnel.pairs,
        hovertemplate="<b>%{y}</b><br>%{x} pairs remain<extra></extra>",
    )
)
fig.update_yaxes(autorange="reversed", ticksuffix="  ")
fig.update_xaxes(title_text="Pairs remaining", showgrid=True, gridcolor=viz.GRID)
viz.style(fig, height=330)
viz.rounded_bars(fig)
viz.value_labels(fig)
st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

st.caption(
    "Criteria are applied in order, so each bar is the survivors of everything "
    "above it. The big fall is the first one: agreement between two independent "
    "labs, in every replicate."
)

# --------------------------------------------------------------------------- #
# near misses
# --------------------------------------------------------------------------- #

st.divider()
st.subheader("The near misses")
st.caption(
    "Pairs that both labs called a hit, every replicate, and still did not "
    "count. Each row shows which later criterion stopped it."
)

near = p[p.c1_hit_both_labs & ~p.passes].copy()


def blocked_by(row: pd.Series) -> str:
    reasons = [label for col, label in data.CRITERIA[1:] if not row[col]]
    return "; ".join(reasons)


if len(near):
    near["Stopped by"] = near.apply(blocked_by, axis=1)
    show = near[
        [
            "name",
            "uuid",
            "toxin",
            "adaptyv_KD_nM",
            "twist_KD_nM",
            "rmax_ratio",
            "Stopped by",
        ]
    ].rename(
        columns={
            "name": "Name",
            "uuid": "ID",
            "toxin": "Toxin",
            "adaptyv_KD_nM": "Adaptyv KD (nM)",
            "twist_KD_nM": "Twist KD (nM)",
            "rmax_ratio": "Rmax ratio",
        }
    )
    st.dataframe(
        show.sort_values("Stopped by"),
        hide_index=True,
        width="stretch",
        column_config={
            "Adaptyv KD (nM)": st.column_config.NumberColumn(format="%.3g"),
            "Twist KD (nM)": st.column_config.NumberColumn(format="%.3g"),
            "Rmax ratio": st.column_config.NumberColumn(format="%.3f"),
        },
    )
    counts = (
        near["Stopped by"].str.split("; ").explode().value_counts().rename("Pairs")
    )
    st.caption(
        "Most common blocker: "
        + ", ".join(f"**{k}** ({v})" for k, v in counts.items())
    )
else:
    st.info("No near misses.")

# --------------------------------------------------------------------------- #
# lab agreement
# --------------------------------------------------------------------------- #

st.divider()
st.subheader("Do the two labs agree?")
st.caption(
    "Each point is one design–toxin pair measured at both labs. The dashed line "
    "is perfect agreement; tighter binding is toward the lower left."
)

both = p.dropna(subset=["adaptyv_KD_nM", "twist_KD_nM"]).copy()
both = both[(both.adaptyv_KD_nM > 0) & (both.twist_KD_nM > 0)]

fig = go.Figure()
groups = [
    ("Other pairs", both[~both.passes], viz.ORANGE),
    ("Binders", both[both.passes], viz.BLUE),
]
for label, sub, color in groups:
    if not len(sub):
        continue
    fig.add_trace(
        go.Scatter(
            x=sub.adaptyv_KD_nM,
            y=sub.twist_KD_nM,
            mode="markers",
            name=label,
            marker=dict(
                size=10, color=color, line=dict(width=2, color=viz.SURFACE), opacity=0.9
            ),
            customdata=sub[["name", "uuid", "toxin"]],
            hovertemplate=(
                "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
                "%{customdata[2]}<br>Adaptyv %{x:.3g} nM · Twist %{y:.3g} nM"
                "<extra></extra>"
            ),
        )
    )

lo = float(min(both.adaptyv_KD_nM.min(), both.twist_KD_nM.min())) * 0.6
hi = float(max(both.adaptyv_KD_nM.max(), both.twist_KD_nM.max())) * 1.6
fig.add_trace(
    go.Scatter(
        x=[lo, hi],
        y=[lo, hi],
        mode="lines",
        line=dict(color=viz.INK_MUTED, width=1, dash="dash"),
        hoverinfo="skip",
        showlegend=False,
    )
)
fig.update_xaxes(type="log", title_text="Adaptyv Bio KD (nM)", showgrid=True, gridcolor=viz.GRID)
fig.update_yaxes(type="log", title_text="Twist Bioscience KD (nM)")
viz.style(fig, height=430, showlegend=True)
st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

if len(both) > 1:
    ratio = both.adaptyv_KD_nM / both.twist_KD_nM
    st.caption(
        f"Across {len(both)} pairs the two labs' KDs differ by a median factor of "
        f"**{ratio.median():.1f}×** (range {ratio.min():.1f}× to {ratio.max():.1f}×). "
        "Adaptyv's values run systematically looser, but the ranking of designs "
        "largely survives the move between labs."
    )

with st.expander("Show every pair and every criterion"):
    st.dataframe(p, hide_index=True, width="stretch", height=400)
    st.download_button(
        "Download as CSV",
        p.to_csv(index=False).encode(),
        file_name="pairs_with_criteria.csv",
        mime="text/csv",
    )
