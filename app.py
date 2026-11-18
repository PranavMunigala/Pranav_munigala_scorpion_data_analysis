"""Entry point: sets up the sidebar navigation and runs the selected page.

The order of PAGES is the order of the sidebar, so the guide sits at the top and
opens by default.
"""

from __future__ import annotations

import streamlit as st

from scorpion import data

st.set_page_config(
    page_title="Scorpion toxin binders",
    page_icon="🦂",
    layout="wide",
    initial_sidebar_state="expanded",
)

PAGES = [
    st.Page(
        "views/start_here.py",
        title="Start here",
        icon="🧭",
        url_path="start",
        default=True,
    ),
    st.Page("views/overview.py", title="Overview", icon="🦂", url_path="overview"),
    st.Page("views/designs.py", title="Designs", icon="🧬", url_path="designs"),
    st.Page(
        "views/binder_criteria.py",
        title="Binder criteria",
        icon="✅",
        url_path="criteria",
    ),
    st.Page(
        "views/in_silico_scores.py",
        title="In silico scores",
        icon="📈",
        url_path="scores",
    ),
    st.Page("views/kinetics.py", title="Kinetics", icon="📊", url_path="kinetics"),
    st.Page(
        "views/visualization.py",
        title="Visualization",
        icon="🔬",
        url_path="structures",
    ),
]

page = st.navigation(PAGES)

st.sidebar.success(f"Data source: {data.source_label()}")
st.sidebar.caption(
    "Dataset: Anthropic/scorpion-toxin-binder-design (CC BY 4.0). "
    "This app computes everything from the released tables."
)

page.run()
