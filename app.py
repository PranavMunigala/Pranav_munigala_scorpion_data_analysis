"""Scorpion toxin binder design — an interactive reading of the data release."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from scorpion import data, viz

st.set_page_config(
    page_title="Scorpion toxin binders",
    page_icon="🦂",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🦂 Scorpion toxin binder design")
st.caption(
    "An interactive reading of Anthropic's data release: 211 proteins designed "
    "*de novo* against five scorpion toxins, screened at two independent labs."
)

d = data.designs()
p = data.pairs()
binders = d[data.is_binder(d.binder_against)]

st.markdown(
    """
Claude designed proteins from scratch and asked whether they would stick to five
scorpion neurotoxins. Two contract research organizations — **Adaptyv Bio** and
**Twist Bioscience** — measured binding independently, without knowing each
other's answers. A design only counts as a *binder* if both labs agree, every
time, and it does not stick to unrelated proteins.

**Seventeen designs cleared that bar.** The pages in the sidebar let you look at
every design, watch the 665 attempted pairings narrow down to those 17, and ask
whether the computational scores predicted any of it.
"""
)

c = st.columns(5)
c[0].metric("Designs", len(d))
c[1].metric("Design–toxin pairs tried", len(p))
c[2].metric("Binders", len(binders))
c[3].metric("Toxins targeted", len(data.TOXINS))
c[4].metric("Independent labs", 2)

st.divider()

# --------------------------------------------------------------------------- #
# binders per toxin
# --------------------------------------------------------------------------- #

left, right = st.columns([3, 2], gap="large")

with left:
    st.subheader("How each toxin fared")
    st.caption(
        "Ts1 defeated every attempt, despite being the second most-targeted toxin."
    )

    per_toxin = (
        p.groupby("toxin")
        .agg(attempted=("uuid", "size"), binders=("published_binder", "sum"))
        .reindex(data.TOXIN_ORDER)
        .reset_index()
    )
    per_toxin["rate"] = 100 * per_toxin.binders / per_toxin.attempted

    fig = go.Figure(
        go.Bar(
            x=per_toxin.toxin,
            y=per_toxin.binders,
            marker_color=viz.BLUE,
            text=per_toxin.binders,
            hovertemplate=(
                "<b>%{x}</b><br>%{y} binders out of %{customdata[0]} designs tried"
                "<br>%{customdata[1]:.1f}% success<extra></extra>"
            ),
            customdata=per_toxin[["attempted", "rate"]],
        )
    )
    fig.update_yaxes(title_text="Binders")
    viz.style(fig, height=320)
    viz.rounded_bars(fig)
    viz.value_labels(fig)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

with right:
    st.subheader("The same numbers, exactly")
    show = per_toxin.rename(
        columns={
            "toxin": "Toxin",
            "attempted": "Designs tried",
            "binders": "Binders",
            "rate": "Success",
        }
    )
    st.dataframe(
        show,
        hide_index=True,
        width="stretch",
        column_config={"Success": st.column_config.NumberColumn(format="%.1f%%")},
    )
    st.caption(
        "A design aimed at several toxins at once appears under each of them, so "
        "the column sums to more than 211."
    )

st.divider()

# --------------------------------------------------------------------------- #
# success by campaign
# --------------------------------------------------------------------------- #

st.subheader("Asking for more at once made success rarer")
st.caption(
    "Later design sets asked a single protein to bind several toxins at once. "
    "The share of designs that produced a binder fell as the ask grew."
)

bd = d.assign(is_binder=data.is_binder(d.binder_against))
bd["group"] = data.design_group(bd.design_set)
grp = (
    bd.groupby("group")
    .agg(designs=("uuid", "size"), binders=("is_binder", "sum"))
    .reindex(data.GROUP_ORDER)
    .reset_index()
)
grp["rate"] = 100 * grp.binders / grp.designs

fig = go.Figure(
    go.Bar(
        x=grp["group"],
        y=grp.rate,
        marker_color=viz.BLUE_RAMP[2],
        text=[f"{r:.0f}%" for r in grp.rate],
        hovertemplate=(
            "<b>%{x}</b><br>%{customdata[0]} of %{customdata[1]} designs "
            "became binders<extra></extra>"
        ),
        customdata=grp[["binders", "designs"]],
    )
)
fig.update_yaxes(title_text="Designs that became binders (%)", ticksuffix="%")
viz.style(fig, height=300)
viz.rounded_bars(fig)
viz.value_labels(fig)
st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

with st.expander("What counts as a binder"):
    st.markdown(
        """
All four of these must hold, for that design against that toxin:

1. **Hit at both labs, in every replicate.** One good run is not enough, and one
   lab is not enough. *(LqhIII was only ever measured at Adaptyv Bio, so it is
   judged on the other three criteria plus a hit there in every replicate.)*
2. **Every main-screen KD at 1 µM or tighter.** KD is the concentration at which
   half the toxin is bound — smaller means stickier.
3. **Fitted Rmax within the theoretical maximum.** A response larger than one
   toxin per binding site can explain means the fit is describing something other
   than clean 1:1 binding.
4. **No CLEC12A signal above background** in either of Twist's two runs. CLEC12A
   is an unrelated human protein: sticking to it means the design is simply
   sticky, not specific.

The **Binder criteria** page re-derives all four from the raw per-replicate data
and checks the result against the published labels.
"""
    )

st.sidebar.success(f"Data source: {data.source_label()}")
st.sidebar.caption(
    "Dataset: Anthropic/scorpion-toxin-binder-design (CC BY 4.0). "
    "This app computes everything from the released tables."
)
