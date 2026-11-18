"""A plain-language guide to what this app shows and how to use it."""

from __future__ import annotations

import streamlit as st

st.set_page_config(page_title="Start here", page_icon="🧭", layout="wide")

st.title("🧭 Start here")
st.caption("What this app is, and how to find your way around it in two minutes.")

st.markdown(
    """
### The one-paragraph version

Claude designed 211 proteins from scratch and asked whether each would stick to a
scorpion neurotoxin. Sticking to a toxin is how an antivenom works: grab the toxin
before it reaches its target in the body. Two laboratories measured every design
independently, without seeing each other's results. **Seventeen designs worked.**
This app lets you look at all of it yourself — not just the successes.
"""
)

st.divider()

# --------------------------------------------------------------------------- #
# the pages
# --------------------------------------------------------------------------- #

st.subheader("The five pages, and what each answers")

pages = [
    (
        "app.py",
        "🦂",
        "Overview",
        "What happened, in numbers",
        "Start here if you want the result and nothing else. How many designs "
        "worked, which toxins proved hard, and why asking one protein to bind "
        "several toxins at once made success rarer.",
    ),
    (
        "pages/1_Designs.py",
        "🧬",
        "Designs",
        "Look up any single design",
        "A searchable list of all 211. Filter by toxin or set, then select a row "
        "to see that design's amino acid sequence, what each lab measured, and "
        "how it was built.",
    ),
    (
        "pages/2_Binder_criteria.py",
        "✅",
        "Binder criteria",
        "Why only 17 survived",
        "A design had to clear four separate bars. This page drops them in one at "
        "a time so you can watch 665 attempts narrow to 17, and shows the near "
        "misses — the ones that failed on a single criterion.",
    ),
    (
        "pages/3_In_silico_scores.py",
        "📈",
        "In silico scores",
        "Could we have predicted it?",
        "Every design was scored by computer before anything was made in a lab. "
        "This page asks how well those predictions matched reality.",
    ),
    (
        "pages/4_Kinetics.py",
        "📊",
        "Kinetics",
        "How the binding actually happens",
        "Binding is a race between how fast a design grabs the toxin and how fast "
        "it lets go. This page plots both for every binder.",
    ),
]

cols = st.columns(len(pages), gap="medium")
for col, (path, icon, name, tagline, body) in zip(cols, pages):
    with col:
        with st.container(border=True):
            st.markdown(f"### {icon}")
            st.markdown(f"**{name}**")
            st.caption(tagline)
            st.write(body)
            st.page_link(path, label=f"Open {name}", icon=icon)

st.divider()

# --------------------------------------------------------------------------- #
# how to drive it
# --------------------------------------------------------------------------- #

left, right = st.columns(2, gap="large")

with left:
    st.subheader("Getting around")
    st.markdown(
        """
- **The sidebar** on the left switches pages. If it is hidden, use the **»**
  arrow at the top-left corner.
- **Charts respond to hover.** Point at any bar or dot to see the exact numbers
  behind it, including the design's name and identifier.
- **Tables are sortable.** Click a column header to sort by it; click again to
  reverse. Drag the edge of a header to widen a column.
- **On the Designs page, click a row** to open that design underneath the table.
- **Anything can be downloaded.** Hover a table and a download icon appears in
  its top-right corner; the Binder criteria page also has an explicit
  *Download as CSV* button for the full per-pair table.
"""
    )

    st.subheader("Three things to try first")
    st.markdown(
        """
1. On **Designs**, turn on *Binders only* and sort by length — the successful
   designs are not all built the same way.
2. On **Binder criteria**, scroll to the scatter plot and find the points far
   from the dashed line. Those are designs the two labs disagreed about most.
3. On **Kinetics**, look at which binders sit furthest toward the upper left.
   Those are the tightest binders in the set.
"""
    )

with right:
    st.subheader("The words you will run into")
    st.markdown(
        """
**Toxin** — the scorpion venom protein we want to neutralize. There are five
here: Cn2, Ts1, AaH II, AaH I and LqhIII.

**Design** — one protein Claude invented, identified by an eight-character code
like `888c8f04`. The ones that worked also got names, like `Cn2-B1`.

**Binder** — a design that genuinely sticks to its toxin, confirmed by both labs
under strict rules. 17 of the 211 qualify.

**Hit** — a single measurement suggesting binding. A hit is not yet a binder: it
has to hold up across replicates and across both labs.

**KD** — binding strength, in nanomolar (nM). It is the concentration at which
half the toxin is caught, so **smaller is stronger**. 1000 nM = 1 µM was the
cutoff for counting at all; the best design here is near 2 nM.

**Replicate** — a repeat of the same measurement. Agreement across replicates is
what separates a real result from a lucky one.

**CLEC12A** — an unrelated human protein used as a trap. A design that sticks to
it is just generally sticky rather than specific to the toxin, so it is
disqualified.

**Rmax** — the largest response physically possible if one toxin binds per site.
A fitted response above it means the curve is describing something other than
clean binding.

**k_on / k_off** — how fast a design attaches to the toxin, and how fast it falls
off. KD is their ratio.

**AUROC** — how well a predicted score sorted the winners from the losers. 0.5 is
a coin flip, 1.0 is perfect.
"""
    )

st.divider()

st.subheader("Two honest caveats")
st.markdown(
    """
- **Small numbers.** 17 successes spread across five toxins means most per-toxin
  and per-set statistics rest on a handful of designs. A single design moving a
  few places can swing a percentage a long way. Treat the AUROC table on the
  *In silico scores* page as suggestive, not settled.
- **Binding is not neutralizing.** Every result here says a design sticks to a
  toxin. Whether that actually protects an animal from venom is a separate
  experiment, and it is not in this dataset.
"""
)

st.info(
    "Everything on these pages is computed live from the released tables — "
    "nothing is hard-coded. The binder criteria are re-derived from the raw "
    "per-replicate measurements and reproduce the published list of 17 exactly.",
    icon="ℹ️",
)
