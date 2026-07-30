#!/usr/bin/env python3
"""Merge classifier, BoltzGen, and Boltz2 scores; rank binders; generate selection plots.

Reads:
  - results/classifiers/merged_classifier_results.csv
  - generations/boltz2_NDM5_KPC3_prediction.csv
  - generations/{KPC3,NDM5}/final_ranked_designs/final_designs_metrics_100.csv

Writes:
  - results/classifiers/merged_binder_selection.csv
  - results/classifiers/selected_binders_top20.csv
  - results/classifiers/plots/*.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parent.parent

BOLTZGEN_METRICS = {
    "KPC3": REPO_ROOT / "generations/KPC3/final_ranked_designs/final_designs_metrics_100.csv",
    "NDM5": REPO_ROOT / "generations/NDM5/final_ranked_designs/final_designs_metrics_100.csv",
}

BG_COLS = [
    "id",
    "final_rank",
    "design_to_target_iptm",
    "design_ptm",
    "min_design_to_target_pae",
    "quality_score",
    "delta_sasa_refolded",
    "plip_hbonds_refolded",
    "filter_rmsd",
    "liability_score",
    "pass_filters",
    "iptm",
    "ptm",
    "min_interaction_pae",
]

B2_COLS = [
    "confidence_score",
    "ptm",
    "iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
]

# Higher composite = better inhibitor-binder candidate (binding up, AMP/toxin/hemo down).
RANK_WEIGHTS = {
    "design_to_target_iptm": 2.0,
    "boltz2_iptm": 2.0,
    "boltz2_confidence": 1.5,
    "quality_score": 1.0,
    "ampscanner_prob_inv": 1.0,
    "macrel_amp_prob_inv": 1.0,
    "macrel_hemo_prob_inv": 0.8,
    "toxin_safe": 1.5,
    "delta_sasa_refolded": 0.5,
}


def _rank_score(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    if series.isna().all():
        return pd.Series(0.0, index=series.index)
    if higher_is_better:
        return series.rank(pct=True, na_option="bottom").fillna(0.0)
    return (1.0 - series.rank(pct=True, na_option="bottom")).fillna(0.0)


def load_boltzgen(target: str) -> pd.DataFrame:
    path = BOLTZGEN_METRICS[target]
    df = pd.read_csv(path)
    keep = [c for c in BG_COLS if c in df.columns]
    out = df[keep + ["sequence"]].copy()
    out = out.rename(columns={c: f"bg_{c}" for c in keep})
    out["target"] = target
    return out


def load_boltz2(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    keep = [c for c in B2_COLS if c in df.columns]
    out = df[["target", "peptide_num"] + keep].copy()
    out = out.rename(columns={c: f"boltz2_{c}" for c in keep})
    return out


def merge_all(
    classifiers_path: Path,
    boltz2_path: Path,
) -> pd.DataFrame:
    clf = pd.read_csv(classifiers_path)
    clf["rank_num"] = clf["seq_id"].str.replace("seq_rank_", "", regex=False).astype(int)

    bg_frames = [load_boltzgen(t) for t in ("KPC3", "NDM5")]
    bg = pd.concat(bg_frames, ignore_index=True)
    bg = bg.rename(columns={"bg_final_rank": "boltzgen_rank"})

    b2 = load_boltz2(boltz2_path)

    df = clf.merge(
        bg,
        on=["target", "sequence"],
        how="left",
        validate="m:1",
    )
    df = df.merge(
        b2,
        left_on=["target", "rank_num"],
        right_on=["target", "peptide_num"],
        how="left",
        suffixes=("", "_dup"),
    )
    if "peptide_num" in df.columns:
        df = df.drop(columns=["peptide_num"])

  # Boolean helpers
    df["macrel_is_amp_bool"] = df["macrel_is_amp"].astype(str).str.lower().eq("true")
    df["is_non_toxin"] = df["toxinpred_class"].eq("Non-Toxin")
    df["is_non_hemo"] = df["macrel_hemolytic"].astype(str).str.lower().ne("hemo")

    return df


def add_composite_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    by_target = out.groupby("target", group_keys=False)

    out["ampscanner_prob_inv"] = 1.0 - out["ampscanner_prob"].fillna(1.0)
    out["macrel_amp_prob_inv"] = 1.0 - out["macrel_amp_prob"].fillna(1.0)
    out["macrel_hemo_prob_inv"] = 1.0 - out["macrel_hemo_prob"].fillna(1.0)
    out["toxin_safe"] = out["is_non_toxin"].astype(float)

    component_cols = []
    for col, higher in [
        ("bg_design_to_target_iptm", True),
        ("boltz2_iptm", True),
        ("boltz2_confidence_score", True),
        ("bg_quality_score", True),
        ("ampscanner_prob_inv", True),
        ("macrel_amp_prob_inv", True),
        ("macrel_hemo_prob_inv", True),
        ("toxin_safe", True),
        ("bg_delta_sasa_refolded", True),
    ]:
        if col not in out.columns:
            continue
        score_col = f"_rank_{col}"
        out[score_col] = by_target[col].transform(lambda s, h=higher: _rank_score(s, h))
        component_cols.append(score_col)

    if component_cols:
        weights = np.array([RANK_WEIGHTS.get(c.replace("_rank_", ""), 1.0) for c in component_cols])
        mat = out[component_cols].to_numpy(dtype=float)
        out["binder_composite_score"] = (mat * weights).sum(axis=1) / weights.sum()
    else:
        out["binder_composite_score"] = 0.0

    out["binder_rank_within_target"] = out.groupby("target")["binder_composite_score"].rank(
        ascending=False, method="first"
    )
    out["safety_score"] = (
        out["toxin_safe"] * 0.4
        + out["macrel_hemo_prob_inv"] * 0.3
        + out["macrel_amp_prob_inv"] * 0.3
    )
    return out


def apply_selection_gates(df: pd.DataFrame) -> pd.Series:
    gates = pd.Series(True, index=df.index)
    gates &= df["is_non_toxin"].fillna(False)
    gates &= df["macrel_hemo_prob"].fillna(1.0) < 0.5
    gates &= ~df["macrel_is_amp_bool"].fillna(False)
    gates &= df["ampscanner_prob"].fillna(1.0) < 0.7
    if "bg_design_to_target_iptm" in df.columns:
        gates &= df["bg_design_to_target_iptm"].fillna(0.0) >= 0.15
    if "boltz2_iptm" in df.columns:
        gates &= df["boltz2_iptm"].fillna(0.0) >= 0.5
    return gates


def select_top(df: pd.DataFrame, per_target: int = 20) -> pd.DataFrame:
    gated = df[apply_selection_gates(df)].copy()
    if gated.empty:
        gated = df.copy()
    return (
        gated.sort_values(["target", "binder_composite_score"], ascending=[True, False])
        .groupby("target", group_keys=False)
        .head(per_target)
        .reset_index(drop=True)
    )


def _savefig(fig: plt.Figure, out_dir: Path, name: str) -> None:
    path = out_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def plot_all(df: pd.DataFrame, selected: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    # 1) Binding: BoltzGen iPTM vs Boltz2 iPTM
    fig, ax = plt.subplots(figsize=(8, 6))
    for target, sub in df.groupby("target"):
        ax.scatter(
            sub["bg_design_to_target_iptm"],
            sub["boltz2_iptm"],
            label=target,
            alpha=0.65,
            s=35,
        )
    ax.set_xlabel("BoltzGen design_to_target_iptm")
    ax.set_ylabel("Boltz2 iPTM")
    ax.set_title("Binding confidence: BoltzGen vs Boltz2")
    ax.legend()
    _savefig(fig, out_dir, "01_boltzgen_vs_boltz2_iptm.png")

    # 2) Composite score distribution by target
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.violinplot(data=df, x="target", y="binder_composite_score", ax=ax, inner="box")
    ax.set_title("Composite binder score distribution")
    _savefig(fig, out_dir, "02_composite_score_violin.png")

    # 3) AMP risk: ampscanner prob vs macrel AMP prob
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        df["ampscanner_prob"],
        df["macrel_amp_prob"],
        c=df["binder_composite_score"],
        cmap="viridis",
        alpha=0.7,
        s=40,
    )
    ax.set_xlabel("AMP Scanner probability")
    ax.set_ylabel("Macrel AMP probability")
    ax.set_title("AMP risk landscape (color = composite score)")
    plt.colorbar(sc, ax=ax, label="composite score")
    _savefig(fig, out_dir, "03_amp_risk_scatter.png")

    # 4) Toxin vs non-toxin counts
    fig, ax = plt.subplots(figsize=(6, 4))
    counts = df.groupby(["target", "toxinpred_class"]).size().unstack(fill_value=0)
    counts.plot(kind="bar", ax=ax, rot=0)
    ax.set_title("ToxinPred class counts")
    ax.set_ylabel("count")
    _savefig(fig, out_dir, "04_toxinpred_counts.png")

    # 5) Hemolysis probability histogram
    fig, ax = plt.subplots(figsize=(8, 5))
    for target, sub in df.groupby("target"):
        ax.hist(sub["macrel_hemo_prob"].dropna(), bins=20, alpha=0.5, label=target)
    ax.axvline(0.5, color="red", ls="--", label="gate 0.5")
    ax.set_xlabel("Macrel hemolysis probability")
    ax.set_ylabel("count")
    ax.set_title("Hemolysis risk")
    ax.legend()
    _savefig(fig, out_dir, "05_hemolysis_histogram.png")

    # 6) Correlation heatmap (key metrics)
    num_cols = [
        "binder_composite_score",
        "bg_design_to_target_iptm",
        "boltz2_iptm",
        "boltz2_confidence_score",
        "bg_quality_score",
        "ampscanner_prob",
        "macrel_amp_prob",
        "macrel_hemo_prob",
        "toxinpred_hybrid_score",
        "bg_delta_sasa_refolded",
        "length",
    ]
    num_cols = [c for c in num_cols if c in df.columns]
    corr = df[num_cols].corr(numeric_only=True)
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", center=0, ax=ax)
    ax.set_title("Metric correlations")
    _savefig(fig, out_dir, "06_correlation_heatmap.png")

    # 7) Top binders bar chart per target
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, (target, sub) in zip(axes, selected.groupby("target")):
        top = sub.head(15).sort_values("binder_composite_score")
        ax.barh(top["seq_id"], top["binder_composite_score"], color="steelblue")
        ax.set_title(f"Top 15 gated binders — {target}")
        ax.set_xlabel("composite score")
    _savefig(fig, out_dir, "07_top_binders_bar.png")

    # 8) Binding vs AMP risk
    fig, ax = plt.subplots(figsize=(8, 6))
    for target, sub in df.groupby("target"):
        ax.scatter(
            sub["bg_design_to_target_iptm"],
            sub["macrel_amp_prob"],
            label=target,
            alpha=0.65,
            s=35,
        )
    ax.set_xlabel("BoltzGen design_to_target_iptm")
    ax.set_ylabel("Macrel AMP probability")
    ax.set_title("Binding vs AMP-like risk")
    ax.legend()
    _savefig(fig, out_dir, "08_binding_vs_amp_risk.png")

    # 9) Boltz2 confidence vs pLDDT
    fig, ax = plt.subplots(figsize=(8, 6))
    for target, sub in df.groupby("target"):
        ax.scatter(
            sub["boltz2_confidence_score"],
            sub["boltz2_complex_plddt"],
            label=target,
            alpha=0.65,
            s=35,
        )
    ax.set_xlabel("Boltz2 confidence")
    ax.set_ylabel("Boltz2 complex pLDDT")
    ax.set_title("Folding quality (Boltz2)")
    ax.legend()
    _savefig(fig, out_dir, "09_boltz2_confidence_vs_plddt.png")

    # 10) Pareto: binding vs safety
    plot_df = df.copy()
    plot_df["safety_score"] = (
        plot_df["toxin_safe"] * 0.4
        + plot_df["macrel_hemo_prob_inv"] * 0.3
        + plot_df["macrel_amp_prob_inv"] * 0.3
    )
    sel = selected.merge(
        plot_df[["seq_id", "target", "safety_score"]],
        on=["seq_id", "target"],
        how="left",
    ).head(30)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        plot_df["bg_design_to_target_iptm"],
        plot_df["safety_score"],
        c=plot_df["binder_composite_score"],
        cmap="plasma",
        alpha=0.75,
        s=45,
    )
    ax.scatter(
        sel["bg_design_to_target_iptm"],
        sel["safety_score"],
        facecolors="none",
        edgecolors="red",
        s=80,
        linewidths=1.5,
        label="top selected",
    )
    ax.set_xlabel("BoltzGen design_to_target_iptm")
    ax.set_ylabel("Safety score (toxin + low hemo + low AMP)")
    ax.set_title("Pareto view: binding vs safety")
    plt.colorbar(sc, ax=ax, label="composite")
    ax.legend()
    _savefig(fig, out_dir, "10_pareto_binding_vs_safety.png")

    # 11) HydrAMP coverage
    fig, ax = plt.subplots(figsize=(6, 4))
    hydramp_counts = df.groupby("target")["hydramp_scored"].sum()
    total = df.groupby("target").size()
    x = hydramp_counts.index
    ax.bar(x, total - hydramp_counts, label="not scored (>25 aa)", color="lightgray")
    ax.bar(x, hydramp_counts, bottom=total - hydramp_counts, label="HydrAMP scored", color="teal")
    ax.set_title("HydrAMP coverage (length ≤ 25 aa)")
    ax.set_ylabel("sequences")
    ax.legend()
    _savefig(fig, out_dir, "11_hydramp_coverage.png")

    # 12) Selected binders: multi-metric parallel bars (top 10 overall)
    top10 = selected.sort_values("binder_composite_score", ascending=False).head(10)
    if not top10.empty:
        metrics = ["bg_design_to_target_iptm", "boltz2_iptm", "macrel_amp_prob_inv", "macrel_hemo_prob_inv"]
        metrics = [m for m in metrics if m in top10.columns]
        norm = top10[["seq_id", "target"] + metrics].copy()
        for m in metrics:
            norm[m] = _rank_score(norm[m], higher_is_better=True)
        fig, ax = plt.subplots(figsize=(12, 6))
        x = np.arange(len(norm))
        width = 0.2
        for i, m in enumerate(metrics):
            ax.bar(x + i * width, norm[m], width=width, label=m.replace("_", " "))
        ax.set_xticks(x + width * (len(metrics) - 1) / 2)
        ax.set_xticklabels(
            [f"{r.seq_id}\n({r.target})" for r in norm.itertuples()],
            rotation=45,
            ha="right",
            fontsize=8,
        )
        ax.set_ylim(0, 1.05)
        ax.set_title("Top selected binders — normalized metric ranks")
        ax.legend(fontsize=8)
        _savefig(fig, out_dir, "12_top_selected_metrics.png")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--classifiers",
        type=Path,
        default=REPO_ROOT / "results/classifiers/merged_classifier_results.csv",
    )
    parser.add_argument(
        "--boltz2",
        type=Path,
        default=REPO_ROOT / "generations/boltz2_NDM5_KPC3_prediction.csv",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=REPO_ROOT / "results/classifiers/merged_binder_selection.csv",
    )
    parser.add_argument(
        "--selected-csv",
        type=Path,
        default=REPO_ROOT / "results/classifiers/selected_binders_top20.csv",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=REPO_ROOT / "results/classifiers/plots",
    )
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    print("Merging classifier + BoltzGen + Boltz2 scores ...")
    df = merge_all(args.classifiers, args.boltz2)
    df = add_composite_scores(df)
    df["passes_selection_gates"] = apply_selection_gates(df)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv} ({len(df)} rows)")

    selected = select_top(df, per_target=args.top_n)
    selected.to_csv(args.selected_csv, index=False)
    print(f"Wrote {args.selected_csv} ({len(selected)} rows)")

    print(f"Generating plots in {args.plots_dir} ...")
    plot_all(df, selected, args.plots_dir)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
