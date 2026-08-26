"""Shared ranking helpers for the AMP + structure panel.

The AMP-filtered pool already passed AMPScanner and Macrel. Ranking therefore
uses Boltz-2 interface confidence, ipSAE, and active-site overlap, with toxin /
hemolysis / developability as soft penalties.
"""
from __future__ import annotations

from difflib import SequenceMatcher

import pandas as pd

TARGETS = ("NDM5", "KPC3")
COLORS = {"NDM5": "#1f77b4", "KPC3": "#ff7f0e"}

# Percentile-rank weights (within target). Higher weight = more influence.
RANK_HIGHER = (
    ("boltz2_iptm", 2.5),
    ("ipSAE_min", 1.5),
    ("epitope_recall", 1.5),
    ("n_catalytic_contacts", 1.0),
    ("interface_precision", 0.8),
    ("pDockQ", 0.6),
    ("ampscanner_prob", 0.4),
    ("macrel_amp_prob", 0.3),
    ("delivery_proxy", 0.3),
)
RANK_LOWER = (
    ("macrel_hemo_prob", 0.8),
    ("aggregation_proxy", 0.4),
    ("n_liabilities", 0.3),
)
TOXIN_SAFE_WEIGHT = 1.0

DEFAULT_MIN_IPTM = 0.5
DEFAULT_N = 25
DEFAULT_MAX_IDENTITY = 0.8

NUMERIC_COLS = [c for c, _ in RANK_HIGHER] + [c for c, _ in RANK_LOWER] + [
    "boltz2_ptm",
    "boltz2_plddt",
    "ipSAE_max",
    "pDockQ2",
    "LIS",
    "length",
    "net_charge_pH7.4",
    "toxinpred_hybrid_score",
]


def find_col(columns, *names: str) -> str | None:
    cols = list(columns)
    for name in names:
        for c in cols:
            if c == name or c.endswith("." + name):
                return c
    return None


def is_true(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    return str(v).strip().lower() in ("true", "1", "yes")


def as_float(v, default=None):
    try:
        if v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def identity(a: str, b: str) -> float:
    return SequenceMatcher(None, a, b).ratio()


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in NUMERIC_COLS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def add_rank_score(df: pd.DataFrame) -> pd.DataFrame:
    """Add rank_score and binder_rank (1 = best) within each target."""
    out = coerce_numeric(df)
    out["toxin_safe"] = out.get("toxinpred_class", pd.Series("", index=out.index)).eq("Non-Toxin").astype(float)
    if "catalytic_ok" in out.columns:
        out["catalytic_ok_bool"] = out["catalytic_ok"].map(is_true)
    else:
        out["catalytic_ok_bool"] = False

    parts = []
    weights = []
    grouped = out.groupby("target", group_keys=False)

    for col, w in RANK_HIGHER:
        if col not in out.columns:
            continue
        score = grouped[col].transform(lambda s: s.rank(pct=True, na_option="keep")).fillna(0.0)
        parts.append(score)
        weights.append(w)
    for col, w in RANK_LOWER:
        if col not in out.columns:
            continue
        score = grouped[col].transform(
            lambda s: s.rank(pct=True, ascending=False, na_option="keep")
        ).fillna(0.0)
        parts.append(score)
        weights.append(w)
    if TOXIN_SAFE_WEIGHT and "toxin_safe" in out.columns:
        parts.append(out["toxin_safe"])
        weights.append(TOXIN_SAFE_WEIGHT)

    if not parts:
        out["rank_score"] = 0.0
    else:
        w = pd.Series(weights, dtype=float)
        mat = pd.concat(parts, axis=1)
        out["rank_score"] = mat.mul(w.values, axis=1).sum(axis=1) / w.sum()

    out["binder_rank"] = out.groupby("target")["rank_score"].rank(ascending=False, method="first")
    return out


def structure_gates(
    df: pd.DataFrame,
    min_iptm: float = DEFAULT_MIN_IPTM,
    require_catalytic: bool = True,
    require_non_toxin: bool = True,
) -> pd.Series:
    gates = pd.Series(True, index=df.index)
    if require_catalytic:
        if "catalytic_ok_bool" in df.columns:
            gates &= df["catalytic_ok_bool"].fillna(False)
        elif "catalytic_ok" in df.columns:
            gates &= df["catalytic_ok"].map(is_true)
    if min_iptm is not None and "boltz2_iptm" in df.columns:
        gates &= df["boltz2_iptm"].fillna(0.0) >= min_iptm
    if require_non_toxin and "toxinpred_class" in df.columns:
        gates &= df["toxinpred_class"].eq("Non-Toxin")
    return gates


def select_diverse(
    df: pd.DataFrame,
    n: int = DEFAULT_N,
    max_identity: float = DEFAULT_MAX_IDENTITY,
) -> pd.DataFrame:
    """Greedy take of top rank_score rows with pairwise sequence identity < max_identity."""
    ordered = df.sort_values(["rank_score", "boltz2_iptm"], ascending=False)
    picked_idx: list = []
    seqs: list[str] = []
    for idx, rec in ordered.iterrows():
        seq = str(rec.get("sequence") or "").upper()
        if not seq:
            continue
        if any(identity(seq, s) >= max_identity for s in seqs):
            continue
        picked_idx.append(idx)
        seqs.append(seq)
        if len(picked_idx) >= n:
            break
    if not picked_idx:
        return df.iloc[0:0].copy()
    return ordered.loc[picked_idx].reset_index(drop=True)


def select_panel(
    df: pd.DataFrame,
    n: int = DEFAULT_N,
    max_identity: float = DEFAULT_MAX_IDENTITY,
    min_iptm: float = DEFAULT_MIN_IPTM,
    require_catalytic: bool = True,
    require_non_toxin: bool = True,
) -> pd.DataFrame:
    """Prefer gated rows; if fewer than n survive diversity, relax iPTM then toxin.

    catalytic_ok is never dropped unless require_catalytic is False.
    """
    scored = add_rank_score(df)
    scored["passes_structure_gates"] = structure_gates(
        scored, min_iptm=min_iptm, require_catalytic=require_catalytic, require_non_toxin=require_non_toxin
    )

    pools = [scored[scored["passes_structure_gates"]]]
    if min_iptm is not None:
        pools.append(
            scored[
                structure_gates(
                    scored, min_iptm=None, require_catalytic=require_catalytic, require_non_toxin=require_non_toxin
                )
            ]
        )
    if require_non_toxin:
        pools.append(
            scored[structure_gates(scored, min_iptm=None, require_catalytic=require_catalytic, require_non_toxin=False)]
        )
    if not require_catalytic:
        pools.append(scored)

    selected_parts = []
    used_ids: set[str] = set()
    remaining = n
    for pool in pools:
        if remaining <= 0:
            break
        pool = pool[~pool["design_id"].astype(str).isin(used_ids)]
        extra = select_diverse(pool, n=remaining, max_identity=max_identity)
        if extra.empty:
            continue
        selected_parts.append(extra)
        used_ids.update(extra["design_id"].astype(str))
        remaining = n - len(used_ids)

    if not selected_parts:
        return scored.iloc[0:0].copy()
    out = pd.concat(selected_parts, ignore_index=True)
    out["panel_rank"] = range(1, len(out) + 1)
    return out
