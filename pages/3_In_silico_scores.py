"""Did the computational scores predict which designs would work?"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from scorpion import data, viz

st.set_page_config(page_title="In silico scores", page_icon="📈", layout="wide")
st.title("📈 Did the scores see it coming?")
st.caption(
    "Every design was scored computationally before anything was made. This page "
    "asks how well those scores separated the eventual binders from the rest."
)

p = data.pairs()
s = data.scores()

M = p.merge(s, on="uuid", how="left", suffixes=("", "_score"))
M["key"] = M.toxin.map(data.TOXIN_KEY)

PREDICTORS = ["ESMFold2", "ESMFold2Fast", "ProtenixV2"]


def per_row(row: pd.Series, metric: str) -> float:
    if metric == "final":
        return row.get(f"{row.key}_final_score", np.nan)
    vals = [row.get(f"{row.key}_{metric}_{pred}", np.nan) for pred in PREDICTORS]
    vals = [v for v in vals if pd.notna(v)]
    return float(np.mean(vals)) if vals else np.nan


METRICS = {
    "Final score": "final",
    "ipSAE (mean of three predictors)": "ipSAE",
    "scDockQ (mean of three predictors)": "scDockQ",
}

choice = st.radio("Score", list(METRICS), horizontal=True)
metric = METRICS[choice]
M["value"] = M.apply(lambda r: per_row(r, metric), axis=1)

usable = M.dropna(subset=["value"])
overall = data.auroc(usable.published_binder, usable.value)

# --------------------------------------------------------------------------- #
# headline
# --------------------------------------------------------------------------- #

m = st.columns(3)
m[0].metric("AUROC, all pairs", f"{overall:.3f}" if not np.isnan(overall) else "—")
m[1].metric("Binders scored", int(usable.published_binder.sum()))
m[2].metric("Non-binders scored", int((~usable.published_binder).sum()))

st.markdown(
    f"""
**AUROC** is the chance that a randomly chosen binder scored higher than a
randomly chosen non-binder. **0.5 is a coin flip; 1.0 is perfect.** At
**{overall:.3f}**, *{choice}* {'carries real signal' if overall > 0.65 else
'separates the two groups only weakly'} across all {len(usable)} scored pairs.
"""
)

st.divider()

# --------------------------------------------------------------------------- #
# distributions
# --------------------------------------------------------------------------- #

left, right = st.columns([3, 2], gap="large")

with left:
    st.subheader("Score distribution")
    st.caption("Each point is one design–toxin pair; the box marks the quartiles.")

    fig = go.Figure()
    series = [
        ("Non-binders", usable[~usable.published_binder], viz.ORANGE),
        ("Binders", usable[usable.published_binder], viz.BLUE),
    ]
    for label, sub, color in series:
        fig.add_trace(
            go.Box(
                x=sub.value,
                name=label,
                marker=dict(color=color, size=7, line=dict(width=1, color=viz.SURFACE)),
                line=dict(color=color, width=2),
                fillcolor="rgba(0,0,0,0)",
                boxpoints="all",
                jitter=0.5,
                pointpos=0,
                hovertemplate="%{x:.3f}<extra>" + label + "</extra>",
            )
        )
    fig.update_xaxes(title_text=choice, showgrid=True, gridcolor=viz.GRID)
    viz.style(fig, height=380, showlegend=False)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

with right:
    st.subheader("By campaign")
    st.caption(
        "Did sets that scored better actually succeed more often? The scores are "
        "compared against the observed success rate of each design group."
    )
    M["group"] = data.design_group(M.design_set)
    d = data.designs().assign(is_binder=data.is_binder(data.designs().binder_against))
    d["group"] = data.design_group(d.design_set)
    obs = d.groupby("group").agg(designs=("uuid", "size"), binders=("is_binder", "sum"))
    obs["Success"] = 100 * obs.binders / obs.designs
    med = M.groupby("group").value.median().rename("Median score")
    joined = obs.join(med).reindex(data.GROUP_ORDER).reset_index()
    st.dataframe(
        joined[["group", "designs", "binders", "Success", "Median score"]].rename(
            columns={"group": "Group", "designs": "Designs", "binders": "Binders"}
        ),
        hide_index=True,
        width="stretch",
        column_config={
            "Success": st.column_config.NumberColumn(format="%.1f%%"),
            "Median score": st.column_config.NumberColumn(format="%.3f"),
        },
    )

st.divider()

# --------------------------------------------------------------------------- #
# within each single-toxin set
# --------------------------------------------------------------------------- #

st.subheader("Within each single-toxin set")
st.caption(
    "The hardest test: among twelve designs aimed at the same toxin, did the "
    "score put the binders on top? Sets with no binder have no AUROC to compute."
)

rows = []
for design_set, toxin in zip(range(1, 6), data.TOXIN_ORDER):
    sub = M[M.design_set == design_set]
    rows.append(
        dict(
            Set=design_set,
            Toxin=toxin,
            Designs=len(sub),
            Binders=int(sub.published_binder.sum()),
            AUROC=data.auroc(sub.published_binder, sub.value),
        )
    )
within = pd.DataFrame(rows)
st.dataframe(
    within,
    hide_index=True,
    width="stretch",
    column_config={"AUROC": st.column_config.NumberColumn(format="%.3f")},
)

st.info(
    "Read these cautiously: with a handful of binders per set, a single design "
    "moving a few ranks swings the number a long way.",
    icon="ℹ️",
)

with st.expander("What these scores are"):
    st.markdown(
        """
- **ipSAE** scores the predicted *interface* between the design and the toxin —
  how confidently the structure predictor places the two chains against each
  other.
- **scDockQ** is a single-chain docking quality score for that same predicted
  complex.
- **Final score** is the combined ranking score used to choose which designs to
  order and test.

Each structural score is predicted by three systems — ESMFold2, ESMFold2Fast and
ProtenixV2 — and averaged here.
"""
    )
