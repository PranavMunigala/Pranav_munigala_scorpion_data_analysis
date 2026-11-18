"""Loading and derivation for the scorpion toxin binder dataset.

The dataset itself is not committed to this repository. This module looks for a
local copy first and otherwise streams the tables straight from Hugging Face, so
the app runs from a fresh clone with nothing downloaded.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

HF_BASE = (
    "https://huggingface.co/datasets/Anthropic/scorpion-toxin-binder-design/"
    "resolve/main/"
)

# (name as it appears in the data, prefix used in design.csv column names)
TOXINS: list[tuple[str, str]] = [
    ("Cn2", "Cn2"),
    ("Ts1", "Ts1"),
    ("AaH II", "AaHII"),
    ("AaH I", "AaHI"),
    ("LqhIII", "LqhIII"),
]
TOXIN_KEY = dict(TOXINS)
TOXIN_ORDER = [t for t, _ in TOXINS]

#: LqhIII was measured only at Adaptyv Bio, so it is judged on criteria 1, 2 and 4.
ADAPTYV_ONLY = {"LqhIII"}

KD_CUTOFF_NM = 1000.0  # criterion 2: 1 uM or tighter


# --------------------------------------------------------------------------- #
# locating the data
# --------------------------------------------------------------------------- #

def local_root() -> Path | None:
    """The local dataset directory, if there is one."""
    env = os.environ.get("SCORPION_DATA_DIR")
    candidates = [Path(env)] if env else []
    here = Path(__file__).resolve().parent.parent
    candidates += [here / "dataset", Path.cwd() / "dataset"]
    for c in candidates:
        if (c / "design.csv").exists():
            return c
    return None


def source_label() -> str:
    root = local_root()
    return str(root) if root else "Hugging Face (streaming)"


def _location(relpath: str) -> str:
    root = local_root()
    return str(root / relpath) if root else HF_BASE + relpath


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #

_UUID_COLS = ("uuid", "uuid_or_control", "parent_uuid", "domainA_uuid", "domainB_uuid")


@st.cache_data(show_spinner="Loading data...")
def table(relpath: str) -> pd.DataFrame:
    """Read any table in the release, keeping identifiers as text."""
    df = pd.read_csv(_location(relpath), dtype={c: str for c in _UUID_COLS})
    return df


def designs() -> pd.DataFrame:
    return table("design.csv")


def scores() -> pd.DataFrame:
    return table("in_silico_scores.csv")


def toxin_table() -> pd.DataFrame:
    return table("toxins.csv")


def sensor_index() -> pd.DataFrame:
    return table("sensorgrams/INDEX.csv")


def origin() -> pd.DataFrame:
    return table("supporting_tables/origin.csv")


def expression() -> pd.DataFrame:
    return table("supporting_tables/expression.csv")


def bvp_elisa() -> pd.DataFrame:
    return table("BVP_ELISA.csv")


def kinetics(uuid: str) -> pd.DataFrame:
    """Per-replicate kinetic fits for one design or control."""
    return table(f"sensorgrams/{uuid}/{uuid}_kinetics.csv")


def is_binder(series: pd.Series) -> pd.Series:
    """True where `binder_against` names a toxin."""
    s = series.astype("string").str.strip()
    return s.notna() & (s != "")


def intended_list(value: object) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip()]


# --------------------------------------------------------------------------- #
# the four binder criteria, re-derived from the released data
# --------------------------------------------------------------------------- #

def _replicate_calls(uuid: str, toxin: str) -> tuple[list[str], list[str]] | None:
    """Per-replicate main-screen calls from the design's kinetics file.

    `sensorgrams/INDEX.csv` carries one summary `result` per row, which for at
    least one design (9a8cdc7a against Cn2) reads `hit` on all three Adaptyv
    replicates while the kinetics file's `Adaptyv_replicate_result` reads
    hit / no hit / no hit. The per-replicate column is the one criterion 1 means,
    so it wins wherever the two disagree. Returns None if the file is missing.
    """
    try:
        k = kinetics(uuid)
    except Exception:
        return None
    k = k[k.toxin == toxin]
    adaptyv, twist = [], []
    for _, row in k.iterrows():
        if row.lab == "Adaptyv Bio":
            call = row.get("Adaptyv_replicate_result")
            adaptyv.append(str(call) if pd.notna(call) else str(row.result))
        elif row.lab == "Twist Bioscience":
            twist.append(str(row.result))
    return adaptyv, twist


@st.cache_data(show_spinner="Re-deriving binder criteria...")
def pairs() -> pd.DataFrame:
    """One row per intended design-toxin pair, with each criterion evaluated.

    A design is a binder of a toxin when it

    1. hits at both labs in every replicate (Adaptyv Bio only, for LqhIII),
    2. has every main-screen KD at 1 uM or tighter,
    3. has a fitted Rmax at Twist within the theoretical Rmax, and
    4. shows no CLEC12A signal above background in either Twist run.

    Criterion 1 is read from `sensorgrams/INDEX.csv`, then re-checked against the
    per-replicate calls in the kinetics files for any pair that would otherwise
    pass -- so the slow per-design files are only opened where the answer hinges
    on them.
    """
    d = designs()
    main = sensor_index()
    main = main[main.assay == "main screen"]

    rows = []
    for _, r in d.iterrows():
        clec = (
            str(r.CLEC12A_Twist_result_run_1) == "no hit"
            and str(r.CLEC12A_Twist_result_run_2) == "no hit"
        )
        for toxin in intended_list(r.intended_toxins):
            k = TOXIN_KEY[toxin]
            reps = main[(main.uuid_or_control == r.uuid) & (main.toxin == toxin)]
            adaptyv = [str(x) for x in reps[reps.lab == "Adaptyv Bio"].result]
            twist = [str(x) for x in reps[reps.lab == "Twist Bioscience"].result]

            akd = r.get(f"{k}_Adaptyv_KD_nM", np.nan)
            tkd = r.get(f"{k}_Twist_KD_nM", np.nan)
            rmax = r.get(f"{k}_Twist_Rmax_ratio", np.nan)
            kds = [x for x in (akd, tkd) if pd.notna(x)]

            def all_hit(calls: list[str]) -> bool:
                return bool(calls) and all(x == "hit" for x in calls)

            def criterion_1(a: list[str], t: list[str]) -> bool:
                return all_hit(a) if toxin in ADAPTYV_ONLY else (all_hit(a) and all_hit(t))

            c1 = criterion_1(adaptyv, twist)
            if c1:  # only then does the finer per-replicate reading matter
                finer = _replicate_calls(r.uuid, toxin)
                if finer is not None:
                    adaptyv, twist = finer
                    c1 = criterion_1(adaptyv, twist)

            c2 = bool(kds) and all(x <= KD_CUTOFF_NM for x in kds)
            c3 = True if toxin in ADAPTYV_ONLY else (pd.notna(rmax) and rmax <= 1.0)
            c4 = clec

            rows.append(
                dict(
                    uuid=r.uuid,
                    name=r["name"],
                    toxin=toxin,
                    design_set=r.design_set,
                    format=r["format"],
                    length_aa=r.length_aa,
                    adaptyv_calls="/".join(adaptyv) or "-",
                    twist_calls="/".join(twist) or "-",
                    n_replicates=len(reps),
                    adaptyv_KD_nM=akd,
                    twist_KD_nM=tkd,
                    rmax_ratio=rmax,
                    c1_hit_both_labs=c1,
                    c2_KD_under_1uM=c2,
                    c3_Rmax_in_range=c3,
                    c4_no_CLEC12A=c4,
                    passes=bool(c1 and c2 and c3 and c4),
                    published_binder=str(r.binder_against).strip() == toxin,
                )
            )

    return pd.DataFrame(rows)


CRITERIA = [
    ("c1_hit_both_labs", "Hit at both labs, every replicate"),
    ("c2_KD_under_1uM", "Every KD at 1 µM or tighter"),
    ("c3_Rmax_in_range", "Fitted Rmax within theoretical Rmax"),
    ("c4_no_CLEC12A", "No CLEC12A off-target signal"),
]


def reproduction_check() -> tuple[int, int, int]:
    """(criteria passes, published binder labels, disagreements)."""
    p = pairs()
    return (
        int(p.passes.sum()),
        int(p.published_binder.sum()),
        int((p.passes != p.published_binder).sum()),
    )


def design_group(design_set: pd.Series) -> pd.Series:
    """The four design campaigns, in the order they were run."""
    return pd.Series(
        np.where(
            design_set <= 5,
            "Single toxin (sets 1-5)",
            np.where(
                design_set == 6,
                "Double, beta (set 6)",
                np.where(design_set == 7, "Triple, alpha (set 7)", "All five (sets 8-10)"),
            ),
        ),
        index=design_set.index,
    )


GROUP_ORDER = [
    "Single toxin (sets 1-5)",
    "Double, beta (set 6)",
    "Triple, alpha (set 7)",
    "All five (sets 8-10)",
]


def auroc(labels, values) -> float:
    """Rank-based AUROC; NaN when either class is empty."""
    y = np.asarray(labels, dtype=bool)
    x = np.asarray(values, dtype=float)
    keep = ~np.isnan(x)
    y, x = y[keep], x[keep]
    pos, neg = int(y.sum()), int((~y).sum())
    if pos == 0 or neg == 0:
        return float("nan")
    ranks = pd.Series(x).rank().to_numpy()
    return float((ranks[y].sum() - pos * (pos + 1) / 2) / (pos * neg))
