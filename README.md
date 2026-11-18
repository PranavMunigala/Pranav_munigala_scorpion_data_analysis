# Scorpion toxin binder design — interactive reanalysis

A Streamlit app and a set of verification scripts for Anthropic's scorpion-toxin
binder design data release. The question behind it: **do the published claims hold
up against the raw released data, and did the computational scores predict any of
it?**

Short answer to the first half — yes, exactly. Applying the four binder criteria to
the per-replicate measurements reproduces the published set of 17 binders with zero
disagreements, across all 665 design–toxin pairs.

## Run it

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run app.py
```

That is all that is needed. The app streams the tables straight from Hugging Face
if no local copy is present, so it works from a fresh clone. To use a local copy
instead, put the release in `dataset/` beside `app.py`, or point
`SCORPION_DATA_DIR` at it:

```bash
hf download Anthropic/scorpion-toxin-binder-design --repo-type dataset --local-dir dataset
```

## What the app shows

| Page | What you can see |
|---|---|
| **Overview** | The headline numbers, how each of the five toxins fared, and how the success rate fell as later design sets asked one protein to bind several toxins at once. |
| **Designs** | All 211 designs, filterable by toxin, set, format and name. Open any one for its sequence, both labs' calls and KDs, off-target result, how it was made, expression, and polyreactivity. |
| **Binder criteria** | The funnel from 665 attempted pairings down to 17 binders, one criterion at a time; the near misses with the reason each was stopped; and a lab-agreement scatter of Adaptyv vs. Twist KDs. |
| **In silico scores** | AUROC of the pre-experiment scores against the eventual outcome, overall, per campaign and within each single-toxin set, with the score distributions. |
| **Kinetics** | k<sub>on</sub> against k<sub>off</sub> for every binder replicate with KD isolines, plus the fitted kinetics for any individual design. |

## The dataset

Not committed here. Get it from
[Anthropic/scorpion-toxin-binder-design](https://huggingface.co/datasets/Anthropic/scorpion-toxin-binder-design)
(CC BY 4.0).

211 proteins designed *de novo* against five scorpion toxins (Cn2, Ts1, AaH II,
AaH I, LqhIII), screened by SPR at two contract research organizations, Adaptyv Bio
and Twist Bioscience. 17 are binders — 7 Cn2, 4 AaH II, 4 AaH I, 2 LqhIII, and none
at all against Ts1.

A design counts as a binder of a toxin when it (1) hits at both labs in every
replicate, (2) has every main-screen KD at 1 µM or tighter, (3) has a fitted Rmax at
Twist within the theoretical Rmax, and (4) shows no CLEC12A signal above background
in either Twist run. LqhIII was measured only at Adaptyv Bio, so it is judged on
criteria 1, 2 and 4.

Report: Claude Science, Nathaniel Horwitz, James van Geelen, Kat Kajderowicz, Amir
Shanehsazzadeh. *De novo design of protein binders against scorpion toxins with
Claude.* October 20, 2026.

## One finding worth knowing

`sensorgrams/INDEX.csv` and the per-design kinetics files do not always agree. For
design `9a8cdc7a` against Cn2, the index reports `hit` on all three Adaptyv
replicates, while that design's kinetics file records `Adaptyv_replicate_result` as
hit / no hit / no hit. Criterion 1 asks for a hit in *every* replicate, so the
per-replicate column is the one that decides it — and it is why `9a8cdc7a` is not a
published binder despite looking like one in the index.

Use the index alone and you get 18 binders, not 17. `scorpion/data.py` reads the
index first and re-checks the kinetics files only for pairs that would otherwise
pass, which both gets the right answer and keeps the app fast when streaming.

## Layout

```
app.py                      Overview page and entry point
pages/                      The other four pages
scorpion/data.py            Loading, and the binder criteria re-derivation
scorpion/viz.py             Shared chart styling and the validated palette
my-analysis/check1-6.py     The original exploratory scripts
my-analysis/pairs_final.csv Per-pair criteria table from check5
```

### The original scripts

`my-analysis/check1.py` through `check6.py` are the exploratory pass that came
first: a census of the release, two successive attempts at re-deriving the binder
criteria, an audit of the CLEC12A off-target calls, the definitive reproduction, and
the in-silico score analysis. They run from the repository root and write to
`analysis/`, so create that directory first. The app supersedes them, but they are
kept because they are where the kinetics-vs-index discrepancy above was found.

## Charts

Colors come from a validated categorical palette; no chart uses more than three
categorical series, every series is directly labelled or has a table view beside it,
and status colors always ship with a word rather than carrying meaning on hue alone.
The app pins a light theme in `.streamlit/config.toml`, which is the surface the
palette was validated against.
