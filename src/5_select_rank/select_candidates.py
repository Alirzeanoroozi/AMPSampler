#!/usr/bin/env python3
"""Select a diverse wet-lab panel from structure-scored AMP binders.

Reads results/structure_manifest_<TARGET>.csv (AMP-filtered designs with Boltz-2,
ipSAE, and active-site overlap). Ranks by interface confidence + epitope
coverage, then keeps the top N sequences with pairwise identity < max-identity.

Default hard gates: catalytic_ok, boltz2_iptm >= 0.5, ToxinPred Non-Toxin.
If that set is too small, iPTM then the toxin gate are relaxed; catalytic_ok
is kept unless --allow-no-catalytic.

Usage (from AMPBinderDesign):
  python src/5_select_rank/select_candidates.py
  python src/5_select_rank/select_candidates.py --targets NDM5 --n 25
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
from ranking import (  # noqa: E402
    DEFAULT_MAX_IDENTITY,
    DEFAULT_MIN_IPTM,
    DEFAULT_N,
    TARGETS,
    add_rank_score,
    select_panel,
    structure_gates,
)

PREFERRED_COLS = [
    "design_id",
    "target",
    "sequence",
    "length",
    "rank_score",
    "binder_rank",
    "panel_rank",
    "passes_structure_gates",
    "boltz2_iptm",
    "boltz2_ptm",
    "boltz2_plddt",
    "ipSAE_min",
    "ipSAE_max",
    "pDockQ",
    "pDockQ2",
    "LIS",
    "epitope_recall",
    "interface_precision",
    "n_catalytic_contacts",
    "catalytic_ok",
    "ampscanner_prob",
    "macrel_amp_prob",
    "macrel_hemolytic",
    "macrel_hemo_prob",
    "hydramp_amp_prob",
    "toxinpred_class",
    "net_charge_pH7.4",
    "delivery_proxy",
    "aggregation_proxy",
    "n_liabilities",
    "liabilities",
]


def load_manifest(results_dir: Path, target: str) -> pd.DataFrame:
    path = results_dir / f"structure_manifest_{target}.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found. Run src/4_structure_prediction/build_manifest.py first."
        )
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"{path} is empty")
    if "target" not in df.columns:
        df["target"] = target
    df["target"] = df["target"].fillna(target)
    return df


def ordered_columns(df: pd.DataFrame) -> list[str]:
    front = [c for c in PREFERRED_COLS if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return front + rest


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    with path.open("w") as fh:
        for r in df.itertuples(index=False):
            did = getattr(r, "design_id", "")
            seq = getattr(r, "sequence", "")
            iptm = getattr(r, "boltz2_iptm", "")
            ipsae = getattr(r, "ipSAE_min", "")
            rec = getattr(r, "epitope_recall", "")
            cat = getattr(r, "n_catalytic_contacts", "")
            score = getattr(r, "rank_score", "")
            fh.write(
                f">{did} iptm={iptm} ipsae={ipsae} epitope_recall={rec} "
                f"n_cat={cat} rank_score={score}\n{seq}\n"
            )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--results-dir", type=Path, default=REPO / "results")
    ap.add_argument("--n", type=int, default=DEFAULT_N, help="panel size per target")
    ap.add_argument("--max-identity", type=float, default=DEFAULT_MAX_IDENTITY)
    ap.add_argument("--min-iptm", type=float, default=DEFAULT_MIN_IPTM)
    ap.add_argument(
        "--allow-no-catalytic",
        action="store_true",
        help="allow designs that do not contact the catalytic core",
    )
    ap.add_argument("--allow-toxin", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)
    all_selected = []

    for target in args.targets:
        pool = load_manifest(args.results_dir, target)
        ranked = add_rank_score(pool)
        ranked["passes_structure_gates"] = structure_gates(
            ranked,
            min_iptm=args.min_iptm,
            require_catalytic=not args.allow_no_catalytic,
            require_non_toxin=not args.allow_toxin,
        )
        selected = select_panel(
            pool,
            n=args.n,
            max_identity=args.max_identity,
            min_iptm=args.min_iptm,
            require_catalytic=not args.allow_no_catalytic,
            require_non_toxin=not args.allow_toxin,
        )
        ranked["selected"] = ranked["design_id"].isin(set(selected["design_id"]))
        ranked_path = args.results_dir / f"ranked_{target}.csv"
        ranked[ordered_columns(ranked)].to_csv(ranked_path, index=False)

        sel_csv = args.results_dir / f"selected_{target}.csv"
        sel_fa = args.results_dir / f"selected_{target}.fasta"
        if selected.empty:
            print(f"[{target}] no designs passed selection")
            continue
        selected[ordered_columns(selected)].to_csv(sel_csv, index=False)
        write_fasta(selected, sel_fa)
        n_gate = int(ranked["passes_structure_gates"].sum()) if "passes_structure_gates" in ranked.columns else 0
        n_sel_gate = int(selected["passes_structure_gates"].sum()) if "passes_structure_gates" in selected.columns else 0
        print(
            f"[{target}] pool={len(ranked)} gated={n_gate} "
            f"selected={len(selected)} (gates_ok={n_sel_gate}) "
            f"id<{args.max_identity} -> {sel_csv.relative_to(REPO)}"
        )
        all_selected.append(selected)

    if all_selected:
        panel = pd.concat(all_selected, ignore_index=True)
        panel_path = args.results_dir / "selected_panel.csv"
        panel[ordered_columns(panel)].to_csv(panel_path, index=False)
        print(f"combined panel {len(panel)} -> {panel_path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
