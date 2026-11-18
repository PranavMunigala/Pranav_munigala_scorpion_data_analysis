# Scorpion toxin binder design — independent reanalysis

Verification and analysis scripts for Anthropic's scorpion-toxin binder design data
release. The goal was to check whether the published claims — in particular the set of
17 binders — can be reproduced from the raw released data, and then to ask what the
in-silico scores actually predicted.

## The data

This repo contains **only my scripts and derived outputs**. The dataset itself is not
redistributed here. Get it from Hugging Face:

```bash
hf download Anthropic/scorpion-toxin-binder-design --repo-type dataset --local-dir dataset
```

Dataset: [Anthropic/scorpion-toxin-binder-design](https://huggingface.co/datasets/Anthropic/scorpion-toxin-binder-design) (CC BY 4.0).
Report: Claude Science, Nathaniel Horwitz, James van Geelen, Kat Kajderowicz, Amir
Shanehsazzadeh. *De novo design of protein binders against scorpion toxins with Claude.*
October 20, 2026.

In brief: 211 de novo designs against five scorpion toxins (Cn2, Ts1, AaH II, AaH I,
LqhIII), screened by SPR at two contract research organizations, Adaptyv Bio and Twist
Bioscience. 17 designs are labelled binders — 7 Cn2, 4 AaH II, 4 AaH I, 2 LqhIII, 0 Ts1.

A design counts as a binder of a toxin when it (1) hits at both labs in every replicate,
(2) has every main-screen KD at 1 µM or tighter, (3) has a fitted Rmax at Twist within the
theoretical Rmax, and (4) shows no CLEC12A signal above background in either Twist run.
LqhIII was measured only at Adaptyv Bio, so it is judged on criteria 1, 2 and 4.

## The scripts

Run from the repository root, with the dataset downloaded to `dataset/`.

| Script | What it does |
|---|---|
| `check1.py` | Census: row counts, design sets and formats, binder tallies against the claimed 17, and the KD table for the binders with each one's Adaptyv/Twist ratio. |
| `check2.py` | First re-derivation of the four binder criteria from the `design.csv` summary columns alone, with a disagreement report against the published `binder_against` label, then an audit of why both-lab hits failed (CLEC12A vs. KD vs. Rmax). |
| `check3.py` | Fixes criterion 1 by taking per-replicate results from `sensorgrams/INDEX.csv` rather than the summary columns. Also inspects two anomalies, `9a8cdc7a` and `cf8652bc`. |
| `check4.py` | Examines the CLEC12A off-target calls, especially the "near background" category, and the `KD_note` field for those two anomalies. |
| `check5.py` | The definitive pass. Reads every `*_kinetics.csv` under `sensorgrams/`, uses `Adaptyv_replicate_result` for per-replicate calls, re-derives all four criteria across the 665 intended design–toxin pairs, and writes `pairs_final.csv`. Then reports success rates by design group. |
| `check6.py` | Analysis rather than verification: AUROC of the in-silico scores (`final_score`, mean ipSAE, mean scDockQ) against the binder outcome, pooled and within sets; median score vs. observed success rate per group; where the binders landed in Claude's own `rank_in_set`; and attribution by `design_tool` and `sequence_design_method`. |

## Output

`pairs_final.csv` — all 665 intended design–toxin pairs with each of the four criteria
evaluated per pair, alongside the published label.

## Note on paths

The scripts write their outputs to `analysis/`, while the committed copy of
`pairs_final.csv` lives in `my-analysis/`. Create `analysis/` before running, or adjust
the `to_csv` paths in `check2.py`, `check5.py` and `check6.py`.
