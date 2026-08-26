#!/usr/bin/env python3
"""Join classifier scores with developability descriptors into one manifest per target.

Writes:
  results/manifest_NDM5.csv
  results/manifest_KPC3.csv

Validates merged_classifier_results.csv so unexpected NaNs / non-numeric values
raise ValueError before a broken table is written. HydrAMP/APEX gaps that are
known (unscored sequences) are kept as empty fields, not as the string 'nan'.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

TARGETS = ("NDM5", "KPC3")

REQUIRED_NO_NAN = [
    "seq_id",
    "target",
    "sequence",
    "length",
    "ampscanner_class",
    "ampscanner_prob",
    "macrel_is_amp",
    "macrel_amp_prob",
    "macrel_hemolytic",
    "macrel_hemo_prob",
    "hydramp_scored",
    "apex_scored",
    "toxinpred_class",
    "toxinpred_hybrid_score",
    "toxinpred_ml_score",
]

REQUIRED_NUMERIC = [
    "length",
    "ampscanner_prob",
    "macrel_amp_prob",
    "macrel_hemo_prob",
    "toxinpred_hybrid_score",
    "toxinpred_ml_score",
]

OPTIONAL_NUMERIC = [
    "hydramp_amp_prob",
    "hydramp_mic_prob",
    "apex_mic_mean_uM",
]

DEV_NUMERIC = [
    "length",
    "net_charge_pH7.4",
    "gravy",
    "hydrophobic_moment",
    "aggregation_proxy",
    "cys_count",
    "delivery_proxy",
    "n_liabilities",
]

ALLOWED_TARGETS = {"NDM5", "KPC3"}
ALLOWED_AMPSCANNER = {"AMP", "Non-AMP"}
ALLOWED_HEMO = {"Hemo", "NonHemo"}
ALLOWED_TOXIN = {"Toxin", "Non-Toxin"}


def _examples(df: pd.DataFrame, mask: pd.Series, col: str = "seq_id", n: int = 5) -> str:
    vals = df.loc[mask, col].astype(str).head(n).tolist()
    return ", ".join(vals) if vals else "(none)"


def _require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def _require_no_nan(df: pd.DataFrame, cols: list[str], label: str) -> None:
    bad = []
    for col in cols:
        mask = df[col].isna() | (df[col].astype(str).str.strip().isin(("", "nan", "NaN", "None", "<NA>")))
        # bool/int columns have no empty-string NaNs; isna is enough
        if df[col].dtype == bool or np.issubdtype(df[col].dtype, np.number):
            mask = df[col].isna()
        n = int(mask.sum())
        if n:
            bad.append(f"{col}: {n} NaN/empty (e.g. {_examples(df, mask)})")
    if bad:
        raise ValueError(f"Unexpected NaN in {label}: " + "; ".join(bad))


def _require_numeric(df: pd.DataFrame, cols: list[str], label: str, allow_na: bool = False) -> None:
    for col in cols:
        coerced = pd.to_numeric(df[col], errors="coerce")
        unparseable = coerced.isna() & df[col].notna() & ~df[col].astype(str).str.strip().isin(("", "nan", "NaN"))
        if int(unparseable.sum()):
            raise ValueError(
                f"Non-numeric values in {label}.{col}: "
                f"{df.loc[unparseable, col].astype(str).head(5).tolist()} "
                f"(ids {_examples(df, unparseable)})"
            )
        inf = coerced.isin([np.inf, -np.inf])
        if int(inf.sum()):
            raise ValueError(f"Inf values in {label}.{col} (ids {_examples(df, inf)})")
        if not allow_na:
            empty = coerced.isna()
            if int(empty.sum()):
                raise ValueError(f"NaN in required numeric {label}.{col} (ids {_examples(df, empty)})")


def _is_true(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def validate_classifiers(df: pd.DataFrame) -> None:
    _require_columns(df, REQUIRED_NO_NAN + OPTIONAL_NUMERIC + ["hydramp_skip_reason"], "classifiers")
    if df.empty:
        raise ValueError("merged_classifier_results.csv is empty")

    _require_no_nan(df, REQUIRED_NO_NAN, "classifiers")
    _require_numeric(df, REQUIRED_NUMERIC, "classifiers", allow_na=False)
    _require_numeric(df, OPTIONAL_NUMERIC, "classifiers", allow_na=True)

    bad_target = ~df["target"].isin(ALLOWED_TARGETS)
    if bad_target.any():
        raise ValueError(f"Unknown target values: {df.loc[bad_target, 'target'].unique().tolist()}")

    if df["seq_id"].duplicated().any():
        dups = df.loc[df["seq_id"].duplicated(), "seq_id"].head(5).tolist()
        raise ValueError(f"Duplicate seq_id in classifiers: {dups}")

    empty_seq = df["sequence"].astype(str).str.strip().eq("")
    if empty_seq.any():
        raise ValueError(f"Empty sequences (ids {_examples(df, empty_seq)})")

    bad_amp = ~df["ampscanner_class"].isin(ALLOWED_AMPSCANNER)
    if bad_amp.any():
        raise ValueError(f"Bad ampscanner_class: {df.loc[bad_amp, 'ampscanner_class'].unique().tolist()}")

    bad_hemo = ~df["macrel_hemolytic"].isin(ALLOWED_HEMO)
    if bad_hemo.any():
        raise ValueError(f"Bad macrel_hemolytic: {df.loc[bad_hemo, 'macrel_hemolytic'].unique().tolist()}")

    bad_tox = ~df["toxinpred_class"].isin(ALLOWED_TOXIN)
    if bad_tox.any():
        raise ValueError(f"Bad toxinpred_class: {df.loc[bad_tox, 'toxinpred_class'].unique().tolist()}")

    scored = _is_true(df["hydramp_scored"])
    hydramp_missing = scored & (df["hydramp_amp_prob"].isna() | df["hydramp_mic_prob"].isna())
    if hydramp_missing.any():
        raise ValueError(
            f"hydramp_scored=True but AMP/MIC probability is NaN "
            f"(n={int(hydramp_missing.sum())}, e.g. {_examples(df, hydramp_missing)})"
        )
    skip_missing = (~scored) & df["hydramp_skip_reason"].isna()
    if skip_missing.any():
        raise ValueError(
            f"hydramp_scored=False but hydramp_skip_reason is NaN "
            f"(n={int(skip_missing.sum())}, e.g. {_examples(df, skip_missing)})"
        )

    apex_n = int(df["apex_mic_mean_uM"].notna().sum())
    if apex_n == 0:
        print("warning: apex_mic_mean_uM is entirely missing (APEX produced no values)")


def validate_developability(df: pd.DataFrame, target: str) -> None:
    need = ["design_id"] + DEV_NUMERIC + ["liabilities"]
    _require_columns(df, need, f"developability_{target}")
    if df.empty:
        raise ValueError(f"developability_{target}.csv is empty")
    _require_no_nan(df, ["design_id"] + DEV_NUMERIC, f"developability_{target}")
    _require_numeric(df, DEV_NUMERIC, f"developability_{target}", allow_na=False)
    if df["design_id"].duplicated().any():
        dups = df.loc[df["design_id"].duplicated(), "design_id"].head(5).tolist()
        raise ValueError(f"Duplicate design_id in developability_{target}: {dups}")


def clean_optional(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["hydramp_skip_reason"] = out["hydramp_skip_reason"].fillna("")
    out["liabilities"] = out.get("liabilities", pd.Series("", index=out.index)).fillna("")
    for col in OPTIONAL_NUMERIC:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def build_target_manifest(clf: pd.DataFrame, dev: pd.DataFrame, target: str) -> pd.DataFrame:
    left = clf[clf["target"] == target].copy()
    if left.empty:
        raise ValueError(f"No classifier rows for target {target}")

    validate_developability(dev, target)

    merged = left.merge(
        dev.drop(columns=["length"], errors="ignore"),
        left_on="seq_id",
        right_on="design_id",
        how="left",
        indicator=True,
        validate="1:1",
    )
    unmatched = merged["_merge"] != "both"
    if unmatched.any():
        raise ValueError(
            f"{target}: {int(unmatched.sum())} designs missing from developability "
            f"(e.g. {_examples(merged, unmatched, 'seq_id')})"
        )
    extra = set(dev["design_id"]) - set(left["seq_id"])
    if extra:
        raise ValueError(
            f"{target}: developability has {len(extra)} ids not in classifiers "
            f"(e.g. {list(extra)[:5]})"
        )

    merged = merged.drop(columns=["_merge", "design_id"])
    merged = merged.rename(columns={"seq_id": "design_id"})
    merged = clean_optional(merged)

    col_order = [
        "design_id",
        "target",
        "sequence",
        "length",
        "ampscanner_class",
        "ampscanner_prob",
        "macrel_is_amp",
        "macrel_amp_prob",
        "macrel_hemolytic",
        "macrel_hemo_prob",
        "hydramp_scored",
        "hydramp_skip_reason",
        "hydramp_amp_prob",
        "hydramp_mic_prob",
        "apex_scored",
        "apex_mic_mean_uM",
        "toxinpred_class",
        "toxinpred_hybrid_score",
        "toxinpred_ml_score",
        "net_charge_pH7.4",
        "gravy",
        "hydrophobic_moment",
        "aggregation_proxy",
        "cys_count",
        "delivery_proxy",
        "n_liabilities",
        "liabilities",
    ]
    return merged[col_order]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--classifiers",
        type=Path,
        default=REPO / "results" / "merged_classifier_results.csv",
        help="merged classifier table (one row per design)",
    )
    ap.add_argument(
        "--dev-dir",
        type=Path,
        default=REPO / "results",
        help="directory containing developability_<TARGET>.csv",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO / "results",
        help="directory for manifest_<TARGET>.csv",
    )
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    args = ap.parse_args()

    if not args.classifiers.is_file():
        fallback = REPO / "results" / "classifiers" / "merged_classifier_results.csv"
        if fallback.is_file():
            args.classifiers = fallback
        else:
            raise FileNotFoundError(args.classifiers)

    clf = pd.read_csv(args.classifiers)
    print(f"Loaded classifiers: {args.classifiers} ({len(clf)} rows)")
    validate_classifiers(clf)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for target in args.targets:
        dev_path = args.dev_dir / f"developability_{target}.csv"
        if not dev_path.is_file():
            raise FileNotFoundError(dev_path)
        dev = pd.read_csv(dev_path)
        manifest = build_target_manifest(clf, dev, target)
        out = args.out_dir / f"manifest_{target}.csv"
        manifest.to_csv(out, index=False, na_rep="")
        leftover_nan = int(manifest.drop(columns=OPTIONAL_NUMERIC, errors="ignore").isna().sum().sum())
        if leftover_nan:
            raise ValueError(f"{out.name} still contains {leftover_nan} unexpected NaN cells")
        print(f"[{target}] {len(manifest)} rows x {len(manifest.columns)} cols -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
