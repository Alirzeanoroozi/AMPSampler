#!/usr/bin/env python3
"""Plots for the structure-scored AMP-binder pool and the selected panel.

Reads results/ranked_<T>.csv (or structure_manifest_<T>.csv) and
results/selected_<T>.csv. Writes PNGs under results/plots/select/.

Usage (from AMPBinderDesign):
  python src/5_select_rank/analyze_binder_selection.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
from ranking import (  # noqa: E402
    COLORS,
    TARGETS,
    add_rank_score,
    is_true,
)

CORR_COLS = [
    "rank_score",
    "boltz2_iptm",
    "ipSAE_min",
    "epitope_recall",
    "interface_precision",
    "n_catalytic_contacts",
    "pDockQ",
    "LIS",
    "ampscanner_prob",
    "macrel_hemo_prob",
    "length",
    "delivery_proxy",
    "aggregation_proxy",
]


def _violin(ax, df: pd.DataFrame, y: str) -> None:
    kw = dict(
        data=df,
        x="target",
        y=y,
        hue="target",
        palette=COLORS,
        inner="quartile",
        cut=0,
        legend=False,
        ax=ax,
    )
    try:
        sns.violinplot(**kw, density_norm="width")
    except TypeError:
        sns.violinplot(**kw, scale="width")


def _style() -> None:
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "savefig.dpi": 200,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 8,
            "legend.frameon": False,
        }
    )


def load_pool(results_dir: Path, targets: list[str]) -> pd.DataFrame:
    frames = []
    for target in targets:
        ranked = results_dir / f"ranked_{target}.csv"
        struct = results_dir / f"structure_manifest_{target}.csv"
        path = ranked if ranked.is_file() else struct
        if not path.is_file():
            raise FileNotFoundError(path)
        df = pd.read_csv(path)
        if "target" not in df.columns:
            df["target"] = target
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    if "rank_score" not in out.columns:
        out = add_rank_score(out)
    if "catalytic_ok_bool" not in out.columns:
        if "catalytic_ok" in out.columns:
            out["catalytic_ok_bool"] = out["catalytic_ok"].map(is_true)
        else:
            out["catalytic_ok_bool"] = False
    if "selected" not in out.columns:
        out["selected"] = False
    else:
        out["selected"] = out["selected"].map(is_true)
    return out


def load_selected(results_dir: Path, targets: list[str]) -> pd.DataFrame:
    frames = []
    for target in targets:
        path = results_dir / f"selected_{target}.csv"
        if path.is_file():
            df = pd.read_csv(path)
            if "target" not in df.columns:
                df["target"] = target
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _savefig(fig: plt.Figure, out_dir: Path, name: str) -> None:
    path = out_dir / name
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(REPO)}")


def _mark_selected(ax, sel: pd.DataFrame, x: str, y: str) -> None:
    if sel is None or sel.empty or x not in sel.columns or y not in sel.columns:
        return
    ax.scatter(
        sel[x],
        sel[y],
        facecolors="none",
        edgecolors="#111111",
        s=90,
        linewidths=1.3,
        zorder=5,
        label="selected",
    )


def plot_iptm_violin(df: pd.DataFrame, sel: pd.DataFrame, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
    _violin(ax, df, "boltz2_iptm")
    if not sel.empty and "boltz2_iptm" in sel.columns:
        sns.stripplot(
            data=sel,
            x="target",
            y="boltz2_iptm",
            color="#111111",
            size=6,
            jitter=0.12,
            ax=ax,
            label="selected",
        )
    ax.axhline(0.5, color="#666666", ls="--", lw=1, label="iPTM = 0.5")
    ax.set_ylabel("Boltz-2 iPTM")
    ax.set_xlabel("")
    ax.set_title("Interface confidence (AMP-filtered pool)")
    ax.set_ylim(0, 1.02)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[-2:], labels[-2:], loc="lower right")
    if own:
        return fig
    return ax


def plot_ipsae_violin(df: pd.DataFrame, sel: pd.DataFrame, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
    _violin(ax, df, "ipSAE_min")
    if not sel.empty and "ipSAE_min" in sel.columns:
        sns.stripplot(
            data=sel,
            x="target",
            y="ipSAE_min",
            color="#111111",
            size=6,
            jitter=0.12,
            ax=ax,
        )
    ax.set_ylabel("ipSAE_min")
    ax.set_xlabel("")
    ax.set_title("ipSAE_min (asymmetric A↔B, cutoffs 10 Å)")
    if own:
        return fig
    return ax


def plot_epitope_vs_iptm(df: pd.DataFrame, sel: pd.DataFrame, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for target, sub in df.groupby("target"):
        ok = sub[sub["catalytic_ok_bool"]]
        miss = sub[~sub["catalytic_ok_bool"]]
        ax.scatter(
            miss["boltz2_iptm"],
            miss["epitope_recall"],
            c=COLORS.get(target, "#333"),
            s=22,
            alpha=0.35,
            marker="x",
            linewidths=0.8,
            label=f"{target} no catalytic contact",
        )
        ax.scatter(
            ok["boltz2_iptm"],
            ok["epitope_recall"],
            c=COLORS.get(target, "#333"),
            s=28,
            alpha=0.7,
            label=f"{target} catalytic_ok",
        )
    _mark_selected(ax, sel, "boltz2_iptm", "epitope_recall")
    ax.axvline(0.5, color="#666666", ls="--", lw=1)
    ax.set_xlabel("Boltz-2 iPTM")
    ax.set_ylabel("Epitope recall")
    ax.set_title("Active-site coverage vs interface confidence")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="best", fontsize=7)
    if own:
        return fig
    return ax


def plot_gate_funnel(df: pd.DataFrame, sel: pd.DataFrame, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
    stages = ["AMP-filtered", "catalytic_ok", "iPTM ≥ 0.5", "Non-Toxin", "selected"]
    x = np.arange(len(stages))
    width = 0.36
    targets = [t for t in TARGETS if t in set(df["target"])]
    for i, target in enumerate(targets):
        sub = df[df["target"] == target]
        n0 = len(sub)
        n1 = int(sub["catalytic_ok_bool"].sum())
        n2 = int((sub["catalytic_ok_bool"] & (sub["boltz2_iptm"] >= 0.5)).sum())
        non_tox = (
            sub["toxinpred_class"].eq("Non-Toxin")
            if "toxinpred_class" in sub.columns
            else pd.Series(True, index=sub.index)
        )
        n3 = int((sub["catalytic_ok_bool"] & (sub["boltz2_iptm"] >= 0.5) & non_tox).sum())
        n4 = 0 if sel.empty else int((sel["target"] == target).sum())
        heights = [n0, n1, n2, n3, n4]
        offset = (i - (len(targets) - 1) / 2) * width
        bars = ax.bar(x + offset, heights, width=width, color=COLORS.get(target), label=target)
        for bar, h in zip(bars, heights):
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1, str(h), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(stages, rotation=15, ha="right")
    ax.set_ylabel("designs")
    ax.set_title("Selection funnel")
    ax.legend(loc="upper right")
    sns.despine(ax=ax)
    if own:
        return fig
    return ax


def plot_hemo_vs_iptm(df: pd.DataFrame, sel: pd.DataFrame, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for target, sub in df.groupby("target"):
        ax.scatter(
            sub["boltz2_iptm"],
            sub["macrel_hemo_prob"],
            c=COLORS.get(target, "#333"),
            s=28,
            alpha=0.65,
            label=target,
        )
    _mark_selected(ax, sel, "boltz2_iptm", "macrel_hemo_prob")
    ax.axhline(0.5, color="#c62828", ls="--", lw=1, label="hemo P = 0.5")
    ax.axvline(0.5, color="#666666", ls="--", lw=1)
    ax.set_xlabel("Boltz-2 iPTM")
    ax.set_ylabel("Macrel hemolysis probability")
    ax.set_title("Binding vs hemolysis risk")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="best", fontsize=8)
    if own:
        return fig
    return ax


def plot_correlation(df: pd.DataFrame, out_dir: Path) -> None:
    cols = [c for c in CORR_COLS if c in df.columns]
    cols = [c for c in cols if pd.to_numeric(df[c], errors="coerce").notna().any()]
    if len(cols) < 2:
        return
    corr = df[cols].corr(numeric_only=True)
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", center=0, ax=ax, square=True, annot_kws={"size": 7})
    ax.set_title("Metric correlations (AMP-filtered pool)")
    _savefig(fig, out_dir, "05_correlation.png")


def plot_catalytic_counts(df: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    counts = (
        df.assign(catalytic=np.where(df["catalytic_ok_bool"], "catalytic_ok", "no catalytic contact"))
        .groupby(["target", "catalytic"])
        .size()
        .unstack(fill_value=0)
    )
    counts.plot(kind="bar", ax=ax, rot=0, color=["#888888", "#2e7d32"])
    ax.set_ylabel("designs")
    ax.set_xlabel("")
    ax.set_title("Catalytic-core contact")
    ax.legend(title="")
    sns.despine(ax=ax)
    _savefig(fig, out_dir, "10_catalytic_ok.png")


def plot_length_vs_iptm(df: pd.DataFrame, sel: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for target, sub in df.groupby("target"):
        ax.scatter(sub["length"], sub["boltz2_iptm"], c=COLORS.get(target), s=28, alpha=0.65, label=target)
    _mark_selected(ax, sel, "length", "boltz2_iptm")
    ax.axhline(0.5, color="#666666", ls="--", lw=1)
    ax.set_xlabel("Binder length (aa)")
    ax.set_ylabel("Boltz-2 iPTM")
    ax.set_title("Length vs interface confidence")
    ax.legend(loc="best")
    _savefig(fig, out_dir, "08_length_vs_iptm.png")


def plot_top_selected(sel: pd.DataFrame, ax=None):
    own = ax is None
    if sel.empty:
        return None if own else ax
    targets = [t for t in TARGETS if t in set(sel["target"])]
    if own:
        fig, axes = plt.subplots(1, len(targets), figsize=(7.0 * len(targets), 6.6), squeeze=False)
        axes = axes[0]
    else:
        fig = None
        axes = [ax]
    for i, target in enumerate(targets):
        a = axes[i] if own else ax
        sub = sel[sel["target"] == target].sort_values("rank_score", ascending=True)
        a.barh(sub["design_id"], sub["rank_score"], color=COLORS.get(target), height=0.7)
        a.set_title(f"Selected panel — {target} (n={len(sub)})")
        a.set_xlabel("composite rank_score")
        a.tick_params(axis="y", labelsize=7)
    if own:
        fig.tight_layout()
        return fig
    return ax


def plot_selected_metrics(sel: pd.DataFrame, out_dir: Path) -> None:
    if sel.empty:
        return
    metrics = [
        ("boltz2_iptm", "iPTM"),
        ("ipSAE_min", "ipSAE_min"),
        ("epitope_recall", "epitope recall"),
        ("interface_precision", "interface precision"),
    ]
    metrics = [(c, lab) for c, lab in metrics if c in sel.columns]
    if not metrics:
        return
    targets = [t for t in TARGETS if t in set(sel["target"])]
    fig, axes = plt.subplots(1, len(targets), figsize=(7.4 * max(len(targets), 1), 7.2), squeeze=False)
    axes = axes[0]
    col_labels = [lab for _, lab in metrics]
    for ax, target in zip(axes, targets):
        sub = sel[sel["target"] == target].sort_values("rank_score", ascending=False)
        mat = sub[[c for c, _ in metrics]].copy()
        mat.index = sub["design_id"].astype(str)
        mat.columns = col_labels
        color = mat.copy()
        for c in color.columns:
            lo, hi = color[c].min(), color[c].max()
            color[c] = (color[c] - lo) / (hi - lo) if hi > lo else 0.5
        sns.heatmap(
            color,
            ax=ax,
            cmap="YlOrRd",
            annot=mat,
            fmt=".2f",
            annot_kws={"size": 7},
            cbar_kws={"shrink": 0.7, "label": "within-panel min–max"},
            vmin=0,
            vmax=1,
        )
        ax.set_title(target)
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=25, labelsize=8)
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle("Selected panel — interface and epitope scores", y=1.01)
    fig.tight_layout()
    _savefig(fig, out_dir, "09_selected_metrics.png")


def plot_overview(df: pd.DataFrame, sel: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.6))
    plot_iptm_violin(df, sel, ax=axes[0, 0])
    plot_ipsae_violin(df, sel, ax=axes[0, 1])
    plot_epitope_vs_iptm(df, sel, ax=axes[0, 2])
    plot_gate_funnel(df, sel, ax=axes[1, 0])
    plot_hemo_vs_iptm(df, sel, ax=axes[1, 1])
    if not sel.empty:
        # compact selected scores: iPTM by target
        sns.stripplot(
            data=sel,
            x="target",
            y="boltz2_iptm",
            hue="target",
            palette=COLORS,
            size=8,
            jitter=0.18,
            legend=False,
            ax=axes[1, 2],
        )
        axes[1, 2].axhline(0.5, color="#666666", ls="--", lw=1)
        axes[1, 2].set_ylim(0, 1.02)
        axes[1, 2].set_xlabel("")
        axes[1, 2].set_ylabel("Boltz-2 iPTM")
        axes[1, 2].set_title("Selected panel iPTM")
    else:
        axes[1, 2].axis("off")
    for ax in axes.ravel():
        sns.despine(ax=ax)
    fig.suptitle("AMP-binder selection after Boltz-2, ipSAE, and active-site overlap", fontsize=14, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _savefig(fig, out_dir, "00_overview.png")


def plot_all(df: pd.DataFrame, sel: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _style()
    plot_overview(df, sel, out_dir)

    fig = plot_iptm_violin(df, sel)
    _savefig(fig, out_dir, "01_iptm_violin.png")
    fig = plot_ipsae_violin(df, sel)
    _savefig(fig, out_dir, "02_ipsae_violin.png")
    fig = plot_epitope_vs_iptm(df, sel)
    _savefig(fig, out_dir, "03_epitope_vs_iptm.png")
    fig = plot_gate_funnel(df, sel)
    _savefig(fig, out_dir, "04_gate_funnel.png")
    plot_correlation(df, out_dir)
    fig = plot_hemo_vs_iptm(df, sel)
    _savefig(fig, out_dir, "06_hemo_vs_iptm.png")
    fig = plot_top_selected(sel)
    if fig is not None:
        _savefig(fig, out_dir, "07_top_selected.png")
    plot_length_vs_iptm(df, sel, out_dir)
    plot_selected_metrics(sel, out_dir)
    plot_catalytic_counts(df, out_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, default=REPO / "results")
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--plots-dir", type=Path, default=REPO / "results" / "plots" / "select")
    args = ap.parse_args()

    df = load_pool(args.results_dir, args.targets)
    sel_ids = set()
    sel = load_selected(args.results_dir, args.targets)
    if not sel.empty:
        sel_ids = set(sel["design_id"].astype(str))
        df["selected"] = df["design_id"].astype(str).isin(sel_ids)
    print(f"pool={len(df)} selected={len(sel)} -> {args.plots_dir}")
    plot_all(df, sel, args.plots_dir)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
