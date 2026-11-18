"""The measured binding kinetics, per design and across the binders."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from scorpion import data, viz

st.title("📊 Binding kinetics")
st.caption(
    "KD says how tightly a design holds a toxin. These are the fitted rates "
    "behind it: how fast it finds the toxin, and how fast it lets go."
)

d = data.designs()
binders = d[data.is_binder(d.binder_against)].copy()

st.markdown(
    """
Binding strength is a race between two rates. **k<sub>on</sub>** is how quickly the
design finds and attaches to the toxin; **k<sub>off</sub>** is how quickly it falls
off again. Their ratio is the KD — so the best binders sit toward the
**upper left**: quick to bind, slow to release.
""",
    unsafe_allow_html=True,
)

# --------------------------------------------------------------------------- #
# all binders at once
# --------------------------------------------------------------------------- #

st.subheader("All seventeen binders")


@st.cache_data(show_spinner="Reading kinetic fits...")
def binder_kinetics() -> pd.DataFrame:
    frames = []
    for _, row in binders.iterrows():
        try:
            k = data.kinetics(row.uuid)
        except Exception:
            continue
        k = k[k.toxin == str(row.binder_against).strip()].copy()
        k["uuid"] = row.uuid
        k["name"] = row["name"]
        frames.append(k)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


K = binder_kinetics()
rates = K.dropna(subset=["kon_per_M_s", "koff_per_s"]) if len(K) else pd.DataFrame()

if len(rates):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=rates.kon_per_M_s,
            y=rates.koff_per_s,
            mode="markers",
            marker=dict(
                size=11,
                color=viz.BLUE,
                line=dict(width=2, color=viz.SURFACE),
                opacity=0.9,
            ),
            customdata=rates[["name", "toxin", "lab", "replicate"]],
            hovertemplate=(
                "<b>%{customdata[0]}</b> vs %{customdata[1]}<br>"
                "%{customdata[2]}, replicate %{customdata[3]}<br>"
                "k_on %{x:.2e} /M/s · k_off %{y:.2e} /s<extra></extra>"
            ),
            showlegend=False,
        )
    )
    # KD isolines: koff = KD * kon
    xs = np.array([rates.kon_per_M_s.min() * 0.5, rates.kon_per_M_s.max() * 2])
    for kd_nm, label in [(1.0, "1 nM"), (100.0, "100 nM"), (1000.0, "1 µM")]:
        kd_m = kd_nm * 1e-9
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=kd_m * xs,
                mode="lines",
                line=dict(color=viz.INK_MUTED, width=1, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        fig.add_annotation(
            x=np.log10(xs[1]),
            y=np.log10(kd_m * xs[1]),
            text=f"KD {label}",
            showarrow=False,
            font=dict(size=11, color=viz.INK_MUTED),
            xanchor="right",
            yanchor="bottom",
        )
    fig.update_xaxes(
        type="log", title_text="k_on  (/M/s) — faster to bind →", showgrid=True,
        gridcolor=viz.GRID,
    )
    fig.update_yaxes(type="log", title_text="k_off  (/s) — ↑ faster to let go")
    viz.style(fig, height=440)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
    st.caption(
        "Each point is one replicate of one binder at one lab. Dotted lines mark "
        "constant KD — points on the same line bind equally tightly by different "
        "combinations of speed and stickiness."
    )
else:
    st.info("No fitted rates available in this copy of the data.")

st.divider()

# --------------------------------------------------------------------------- #
# one design
# --------------------------------------------------------------------------- #

st.subheader("One design in detail")

labelled = binders.assign(label=lambda x: x["name"] + "  (" + x.uuid + ")")
others = d[~data.is_binder(d.binder_against)].assign(
    label=lambda x: x.uuid + "  — not a binder"
)
options = labelled.label.tolist() + others.label.tolist()
pick = st.selectbox("Design", options, index=0)
uuid = pick.split("(")[-1].rstrip(")").strip() if "(" in pick else pick.split()[0]

try:
    k = data.kinetics(uuid)
except Exception:
    st.warning(f"No kinetics file for {uuid} in this copy of the data.")
    st.stop()

assay_cols = [
    c
    for c in [
        "lab",
        "toxin",
        "replicate",
        "result",
        "KD_nM",
        "kon_per_M_s",
        "koff_per_s",
        "fit_model",
        "Adaptyv_fit_quality",
        "Adaptyv_replicate_result",
        "method",
    ]
    if c in k.columns
]

st.dataframe(
    k[assay_cols].rename(
        columns={
            "lab": "Lab",
            "toxin": "Toxin",
            "replicate": "Replicate",
            "result": "Result",
            "KD_nM": "KD (nM)",
            "kon_per_M_s": "k_on (/M/s)",
            "koff_per_s": "k_off (/s)",
            "fit_model": "Fit model",
            "Adaptyv_fit_quality": "Fit quality",
            "Adaptyv_replicate_result": "Per-replicate call",
            "method": "Method",
        }
    ),
    hide_index=True,
    width="stretch",
    column_config={
        "KD (nM)": st.column_config.NumberColumn(format="%.3g"),
        "k_on (/M/s)": st.column_config.NumberColumn(format="%.2e"),
        "k_off (/s)": st.column_config.NumberColumn(format="%.2e"),
    },
)

kd = k.dropna(subset=["KD_nM"]).copy()
if len(kd):
    kd["who"] = kd.lab + " · " + kd.toxin + " rep " + kd.replicate.astype(str)
    fig = go.Figure(
        go.Bar(
            x=kd.KD_nM,
            y=kd.who,
            orientation="h",
            marker_color=viz.BLUE,
            text=[f"{v:.3g}" for v in kd.KD_nM],
            hovertemplate="%{y}<br>KD %{x:.3g} nM<extra></extra>",
        )
    )
    fig.update_xaxes(
        type="log", title_text="KD (nM) — tighter binding to the left", showgrid=True,
        gridcolor=viz.GRID,
    )
    fig.update_yaxes(autorange="reversed")
    viz.style(fig, height=60 + 34 * len(kd))
    viz.rounded_bars(fig)
    viz.value_labels(fig)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
else:
    st.caption("No fitted KD for this design — it was not a hit, or the fit failed.")

st.info(
    "This page reads the fitted kinetics tables. The raw sensorgram traces and "
    "images live in the same per-design folders on Hugging Face but are not part "
    "of a CSV-only download.",
    icon="ℹ️",
)
