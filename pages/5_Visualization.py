"""Three-dimensional views of the designs and the toxins they were built to grip."""

from __future__ import annotations

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from scorpion import data, viz

st.set_page_config(page_title="Visualization", page_icon="🔬", layout="wide")
st.title("🔬 Structure visualization")
st.caption(
    "The shapes behind the numbers: how each design was meant to wrap around its "
    "toxin, and where the two touch."
)

d = data.designs()
idx = data.structure_index()
faces = data.predicted_interfaces()

# --------------------------------------------------------------------------- #
# choose a structure
# --------------------------------------------------------------------------- #

binders = d[data.is_binder(d.binder_against)]
binder_ids = binders.uuid.tolist()

c1, c2, c3 = st.columns([3, 2, 3])

with c1:
    only_binders = st.toggle("Binders only", value=True)
    pool = binders if only_binders else d
    labels = {}
    for _, row in pool.iterrows():
        name = row["name"] if pd.notna(row["name"]) else row.uuid
        tag = f" — binds {row.binder_against}" if pd.notna(row.binder_against) else ""
        labels[f"{name} ({row.uuid}){tag}"] = row.uuid
    if not labels:
        st.warning("Nothing to show.", icon="🔍")
        st.stop()
    pick = st.selectbox("Design", list(labels))
    uuid = labels[pick]

mine = idx[idx.uuid == uuid]

with c2:
    kinds = [k for k in ["design", "predicted complex", "predicted monomer"] if k in set(mine.kind)]
    kind = st.selectbox(
        "Structure",
        kinds,
        help=(
            "**design** is the model Claude built. **predicted complex** is an "
            "independent structure predictor's attempt at the same pairing. "
            "**predicted monomer** is the design alone, with no toxin."
        ),
    )
    sub = mine[mine.kind == kind]

with c3:
    cc1, cc2 = st.columns(2)
    toxins_here = sorted({t for t in sub.toxin.dropna().unique()})
    toxin = cc1.selectbox("Toxin", toxins_here) if toxins_here else None
    if toxin:
        sub = sub[sub.toxin == toxin]
    preds = sorted({p for p in sub.predictor.dropna().unique()})
    predictor = cc2.selectbox("Predictor", preds) if preds else None
    if predictor:
        sub = sub[sub.predictor == predictor]

if not len(sub):
    st.warning("No structure matches that combination.", icon="🔍")
    st.stop()

record = sub.iloc[0]

# --------------------------------------------------------------------------- #
# view options
# --------------------------------------------------------------------------- #

o1, o2, o3, o4 = st.columns([2, 2, 2, 2])
coloring = o1.radio(
    "Color by",
    ["Chain", "Confidence"],
    horizontal=True,
    help=(
        "**Chain** separates the design from the toxin. **Confidence** shades each "
        "residue by the predictor's own certainty, stored in the B-factor column — "
        "dark blue is confident, orange is not."
    ),
)
representation = o2.radio("Style", ["Cartoon", "Cartoon + surface"], horizontal=True)
show_contacts = o3.toggle(
    "Highlight contacts",
    value=True,
    help="Draws the residues the two chains touch as sticks.",
)
spin = o4.toggle("Spin", value=False)

try:
    pdb = data.structure_text(record.file)
except Exception as exc:  # network failure, missing file
    st.error(
        f"Could not load `{record.file}`. Structures are not part of a CSV-only "
        f"download, so this page needs either the full dataset locally or an "
        f"internet connection.\n\n`{exc}`",
        icon="⚠️",
    )
    st.stop()

# contact residues, where the release lists them
contacts: dict[str, list[int]] = {"A": [], "B": []}
face = pd.DataFrame()
if toxin:
    face = faces[(faces.uuid == uuid) & (faces.toxin == toxin)]
    if predictor:
        want = predictor.lower().replace(" ", "-").replace("v2", "v2")
        narrowed = face[face.source.str.lower() == want]
        if len(narrowed):
            face = narrowed
if len(face):
    row = face.iloc[0]
    for chain, column in (("A", "design_residues_in_contact"), ("B", "toxin_residues_in_contact")):
        raw = row.get(column)
        if pd.notna(raw):
            contacts[chain] = [int(x) for x in str(raw).split(";") if x.strip().isdigit()]

# --------------------------------------------------------------------------- #
# render
# --------------------------------------------------------------------------- #

view = py3Dmol.view(width=1100, height=560)
view.addModel(pdb, "pdb")
view.setBackgroundColor(viz.SURFACE)

if coloring == "Chain":
    view.setStyle({"chain": "A"}, {"cartoon": {"color": viz.BLUE}})
    view.setStyle({"chain": "B"}, {"cartoon": {"color": viz.ORANGE}})
else:
    view.setStyle(
        {},
        {
            "cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "roygb",
                    "min": 50,
                    "max": 95,
                }
            }
        },
    )

if show_contacts:
    for chain, residues in contacts.items():
        if residues:
            view.addStyle(
                {"chain": chain, "resi": residues},
                {"stick": {"color": viz.AQUA, "radius": 0.18}},
            )

if representation == "Cartoon + surface":
    view.addSurface(
        "VDW",
        {"opacity": 0.55, "color": viz.ORANGE},
        {"chain": "B"},
    )

view.zoomTo()
if spin:
    view.spin(True)

components.html(view._make_html(), height=580, scrolling=False)

# legend — named, never color alone
chains = sorted({line[21] for line in pdb.splitlines() if line.startswith("ATOM")})
swatch = (
    "<span style='display:inline-block;width:11px;height:11px;border-radius:3px;"
    "background:{c};margin-right:6px;vertical-align:middle'></span>"
)
if coloring == "Chain":
    legend = f"{swatch.format(c=viz.BLUE)}<b>Design</b> (chain A)"
    if "B" in chains:
        legend += f" &nbsp;&nbsp; {swatch.format(c=viz.ORANGE)}<b>{toxin or 'Toxin'}</b> (chain B)"
else:
    legend = "Shaded by predictor confidence: <b>blue, confident → orange, uncertain</b>"
if show_contacts and any(contacts.values()):
    legend += f" &nbsp;&nbsp; {swatch.format(c=viz.AQUA)}<b>Contact residues</b> (sticks)"
st.markdown(legend, unsafe_allow_html=True)

st.caption(
    "Drag to rotate, scroll to zoom, right-drag to pan. "
    f"Showing `{record.file}`."
)

st.divider()

# --------------------------------------------------------------------------- #
# what the interface looks like, in numbers
# --------------------------------------------------------------------------- #

left, right = st.columns([2, 3], gap="large")

with left:
    st.subheader("The grip")
    if len(face):
        row = face.iloc[0]
        st.caption(
            "How much surface the two proteins bury against each other. More "
            "buried area generally means a firmer hold."
        )
        m = st.columns(2)
        m[0].metric("Buried area", f"{float(row.buried_area_A2):,.0f} Å²")
        m[1].metric("Apolar share", f"{100 * float(row.apolar_share):.0f}%")
        m = st.columns(2)
        m[0].metric("Design residues touching", int(row.n_design_residues_in_contact))
        m[1].metric("Toxin residues touching", int(row.n_toxin_residues_in_contact))
        if pd.notna(row.hotspots_in_contact):
            st.markdown("**Hotspots gripped**")
            st.code(str(row.hotspots_in_contact), language=None)
            st.caption(
                "Hotspots are the toxin residues known to matter most. A design "
                "that covers them is covering the business end of the toxin."
            )
    else:
        st.info(
            "No interface measurements published for this combination.", icon="ℹ️"
        )

with right:
    st.subheader("Every structure released for this design")
    st.caption(
        f"{len(mine)} files: the design model, plus three predictors' attempts at "
        "each toxin pairing and at the design on its own."
    )
    st.dataframe(
        mine[["kind", "predictor", "toxin", "file"]].rename(
            columns={
                "kind": "Kind",
                "predictor": "Predictor",
                "toxin": "Toxin",
                "file": "File",
            }
        ),
        hide_index=True,
        width="stretch",
        height=300,
    )

st.info(
    "Structure files are not part of a CSV-only download. This page pulls them "
    "from Hugging Face on demand unless a full local copy is present.",
    icon="ℹ️",
)
