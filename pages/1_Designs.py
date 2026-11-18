"""Browse all 211 designs and open any one of them."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from scorpion import data, viz

st.set_page_config(page_title="Designs", page_icon="🧬", layout="wide")
st.title("🧬 Designs")
st.caption("Every design in the release, with its result at each lab.")

d = data.designs()
p = data.pairs()
org = data.origin()

# --------------------------------------------------------------------------- #
# filters
# --------------------------------------------------------------------------- #

f1, f2, f3, f4 = st.columns([2, 2, 2, 3])
toxin_pick = f1.multiselect("Intended toxin", data.TOXIN_ORDER, default=[])
set_pick = f2.multiselect("Design set", sorted(d.design_set.unique()), default=[])
fmt_pick = f3.multiselect("Format", sorted(d["format"].dropna().unique()), default=[])
only_binders = f4.toggle("Binders only", value=False)
query = f4.text_input("Search name or id", placeholder="Cn2-B1, 888c8f04 ...")

view = d.copy()
view["is_binder"] = data.is_binder(view.binder_against)
if toxin_pick:
    view = view[
        view.intended_toxins.apply(
            lambda v: any(t in data.intended_list(v) for t in toxin_pick)
        )
    ]
if set_pick:
    view = view[view.design_set.isin(set_pick)]
if fmt_pick:
    view = view[view["format"].isin(fmt_pick)]
if only_binders:
    view = view[view.is_binder]
if query:
    q = query.strip().lower()
    view = view[
        view.uuid.str.lower().str.contains(q, na=False)
        | view["name"].astype("string").str.lower().str.contains(q, na=False)
    ]

st.caption(f"**{len(view)}** of {len(d)} designs shown. Select a row to open it.")

table = view[
    [
        "uuid",
        "name",
        "design_set",
        "intended_toxins",
        "format",
        "length_aa",
        "binder_against",
    ]
].rename(
    columns={
        "uuid": "ID",
        "name": "Name",
        "design_set": "Set",
        "intended_toxins": "Designed against",
        "format": "Format",
        "length_aa": "Length (aa)",
        "binder_against": "Binds",
    }
)

event = st.dataframe(
    table,
    hide_index=True,
    width="stretch",
    height=380,
    on_select="rerun",
    selection_mode="single-row",
)

if not len(view):
    st.warning("No designs match these filters.", icon="○")
    st.stop()

rows = event.selection.rows if event and event.selection else []
if not rows:
    st.caption("Showing the first design in the list — select any row to switch.")
chosen = view.iloc[rows[0] if rows else 0]

# --------------------------------------------------------------------------- #
# one design in full
# --------------------------------------------------------------------------- #

st.divider()
title = chosen["name"] if pd.notna(chosen["name"]) else chosen.uuid
st.subheader(f"{title}")

if pd.notna(chosen.binder_against):
    st.success(f"Binder of **{chosen.binder_against}**", icon="✅")
else:
    st.warning("Not a binder of any toxin it was designed against", icon="○")

m = st.columns(5)
m[0].metric("ID", chosen.uuid)
m[1].metric("Design set", int(chosen.design_set))
m[2].metric("Length", f"{int(chosen.length_aa)} aa")
m[3].metric("Format", chosen["format"])
m[4].metric("Designed against", chosen.intended_toxins)

st.markdown("**Sequence**")
st.code(str(chosen.sequence), language=None, wrap_lines=True)

col_a, col_b = st.columns([3, 2], gap="large")

with col_a:
    st.markdown("**Measured against each toxin it was designed for**")
    sub = p[p.uuid == chosen.uuid]
    out = sub[
        [
            "toxin",
            "adaptyv_calls",
            "twist_calls",
            "adaptyv_KD_nM",
            "twist_KD_nM",
            "rmax_ratio",
            "passes",
        ]
    ].rename(
        columns={
            "toxin": "Toxin",
            "adaptyv_calls": "Adaptyv replicates",
            "twist_calls": "Twist replicates",
            "adaptyv_KD_nM": "Adaptyv KD (nM)",
            "twist_KD_nM": "Twist KD (nM)",
            "rmax_ratio": "Rmax ratio",
            "passes": "Binder",
        }
    )
    st.dataframe(
        out,
        hide_index=True,
        width="stretch",
        column_config={
            "Adaptyv KD (nM)": st.column_config.NumberColumn(format="%.3g"),
            "Twist KD (nM)": st.column_config.NumberColumn(format="%.3g"),
            "Rmax ratio": st.column_config.NumberColumn(format="%.3f"),
            "Binder": st.column_config.CheckboxColumn(),
        },
    )
    if pd.notna(chosen.KD_note):
        st.caption(f"Note from the release: {chosen.KD_note}")

    st.markdown("**Off-target check — CLEC12A at Twist Bioscience**")
    o1, o2 = st.columns(2)
    o1.metric("Run 1", str(chosen.CLEC12A_Twist_result_run_1))
    o2.metric("Run 2", str(chosen.CLEC12A_Twist_result_run_2))

with col_b:
    st.markdown("**How it was made**")
    o = org[org.uuid == chosen.uuid]
    if len(o):
        o = o.iloc[0]
        made = {
            "Rank in its set": o.rank_in_set,
            "Design tool": o.design_tool,
            "Backbone family": o.backbone_family,
            "Backbone made against": o.backbone_made_against,
            "Sequence design": o.sequence_design_method,
            "Redesign round": o.redesign_round,
        }
        st.dataframe(
            pd.DataFrame(
                {"": list(made), "Value": [str(v) for v in made.values()]}
            ).set_index(""),
            width="stretch",
        )
        st.caption(
            "*Rank in its set* is where Claude itself placed this design before "
            "any experiment was run — 1 was its top pick."
        )

    exp = data.expression()
    e = exp[exp.uuid == chosen.uuid]
    if len(e):
        e = e.iloc[0]
        st.markdown("**Expression**")
        st.metric("Purified yield", f"{e.purified_yield_mg_per_mL} mg/mL")

    bvp = data.bvp_elisa()
    b = bvp[bvp.uuid == chosen.uuid]
    if len(b):
        b = b.iloc[0]
        st.markdown("**Polyreactivity (BVP ELISA)**")
        st.metric(
            "OD450, mean", f"{b.OD450_mean}", help="Lower means less generally sticky."
        )
        st.caption(f"{b.percent_of_cixutumumab}% of the cixutumumab reference.")
