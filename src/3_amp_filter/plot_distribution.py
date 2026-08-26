#!/usr/bin/env python3
"""Bar / histogram distributions of every manifest metric in one PNG.

Each panel has its own axes. Categorical columns are count bars; numeric
columns are histograms. NDM5 and KPC3 are overlaid on every panel.

Usage:
  python src/3_amp_filter/plot_distribution.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

TARGETS = ("NDM5", "KPC3")
COLORS = {"NDM5": "#1f77b4", "KPC3": "#ff7f0e"}

CATEGORICAL = [
    "ampscanner_class",
    "macrel_is_amp",
    "macrel_hemolytic",
    "toxinpred_class",
    "hydramp_scored",
    "hydramp_skip_reason",
]

NUMERIC = [
    "length",
    "ampscanner_prob",
    "macrel_amp_prob",
    "macrel_hemo_prob",
    "hydramp_amp_prob",
    "hydramp_mic_prob",
    "toxinpred_hybrid_score",
    "toxinpred_ml_score",
    "net_charge_pH7.4",
    "gravy",
    "hydrophobic_moment",
    "aggregation_proxy",
    "cys_count",
    "delivery_proxy",
    "n_liabilities",
]


def load_manifests(manifest_dir: Path, targets: list[str]) -> pd.DataFrame:
    frames = []
    for target in targets:
        path = manifest_dir / f"manifest_{target}.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        df = pd.read_csv(path)
        if df.empty:
            raise ValueError(f"{path.name} is empty")
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _plot_categorical(ax, df: pd.DataFrame, col: str) -> None:
    series = df[col].fillna("").astype(str)
    if col == "hydramp_skip_reason":
        series = series.replace({"": "scored"})
    else:
        series = series.replace({"": "(missing)"})
    plot_df = df.copy()
    plot_df[col] = series
    cats = list(dict.fromkeys(series.tolist()))
    x = np.arange(len(cats))
    width = 0.38 if df["target"].nunique() > 1 else 0.7
    n_t = df["target"].nunique()
    for i, (target, sub) in enumerate(plot_df.groupby("target")):
        counts = sub[col].value_counts()
        heights = [int(counts.get(c, 0)) for c in cats]
        offset = (i - (n_t - 1) / 2) * width
        ax.bar(x + offset, heights, width=width, label=target, color=COLORS.get(target), alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, rotation=25, ha="right", fontsize=7)
    ax.set_ylabel("count")
    ax.set_title(col, fontsize=9)
    ax.legend(fontsize=7, loc="best", frameon=False)


def _plot_numeric(ax, df: pd.DataFrame, col: str) -> None:
    vals = []
    for target, sub in df.groupby("target"):
        v = pd.to_numeric(sub[col], errors="coerce").dropna()
        if v.empty:
            continue
        vals.append(v.to_numpy())
        ax.hist(
            v,
            bins=min(30, max(8, int(np.sqrt(len(v))))),
            alpha=0.55,
            label=f"{target} (n={len(v)})",
            color=COLORS.get(target),
            edgecolor="none",
        )
    ax.set_ylabel("count")
    ax.set_title(col, fontsize=9)
    if vals:
        ax.legend(fontsize=7, loc="best", frameon=False)
    else:
        ax.text(0.5, 0.5, "no values", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])


def plot_all(df: pd.DataFrame, out_path: Path) -> None:
    cat_cols = [c for c in CATEGORICAL if c in df.columns]
    num_cols = [c for c in NUMERIC if c in df.columns]
    # drop numeric columns that are entirely missing
    num_cols = [c for c in num_cols if pd.to_numeric(df[c], errors="coerce").notna().any()]
    panels = [("cat", c) for c in cat_cols] + [("num", c) for c in num_cols]
    if not panels:
        raise ValueError("no plottable columns in manifests")

    n = len(panels)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 3.1 * nrows),
        squeeze=False,
        sharex=False,
        sharey=False,
    )
    fig.suptitle("AMP filter metric distributions", fontsize=14, y=0.995)

    for i, (kind, col) in enumerate(panels):
        ax = axes[i // ncols][i % ncols]
        if kind == "cat":
            _plot_categorical(ax, df, col)
        else:
            _plot_numeric(ax, df, col)

    for j in range(len(panels), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path} ({n} panels, {len(df)} designs)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-dir", type=Path, default=REPO / "results")
    ap.add_argument(
        "--out",
        type=Path,
        default=REPO / "results" / "plots" / "amp_filter_distributions.png",
    )
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    args = ap.parse_args()

    df = load_manifests(args.manifest_dir, args.targets)
    plot_all(df, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
