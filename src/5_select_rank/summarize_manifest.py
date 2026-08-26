#!/usr/bin/env python3
"""Quality summary of a structure-scored AMP-binder table.

Default input is results/ranked_<T>.csv (or structure_manifest_<T>.csv).
Reports catalytic-core yield, Boltz-2 / ipSAE / epitope distributions, and
the top-N by rank_score.

Usage (from AMPBinderDesign):
  python src/5_select_rank/summarize_manifest.py
  python src/5_select_rank/summarize_manifest.py --manifest results/ranked_NDM5.csv
"""
from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
from ranking import TARGETS, find_col, is_true  # noqa: E402

KEY_METRICS = [
    ("boltz2_iptm", "Boltz-2 iPTM"),
    ("ipSAE_min", "ipSAE_min (A↔B)"),
    ("epitope_recall", "Active-site coverage"),
    ("interface_precision", "Interface focus on epitope"),
    ("n_catalytic_contacts", "# catalytic-core contacts"),
    ("pDockQ", "pDockQ"),
    ("LIS", "LIS"),
    ("boltz2_plddt", "Boltz-2 complex pLDDT"),
    ("ampscanner_prob", "AMPScanner P(AMP)"),
    ("macrel_amp_prob", "Macrel P(AMP)"),
    ("macrel_hemo_prob", "Macrel P(hemolytic)"),
    ("length", "Binder length (aa)"),
    ("net_charge_pH7.4", "Net charge at pH 7.4"),
    ("delivery_proxy", "Periplasmic-delivery proxy"),
    ("aggregation_proxy", "Aggregation proxy"),
    ("n_liabilities", "# synthesis liabilities"),
    ("rank_score", "Composite rank score"),
]


def numeric_stats(vals):
    nums = []
    for v in vals:
        try:
            if v is None or str(v).strip() == "":
                continue
            nums.append(float(v))
        except (TypeError, ValueError):
            pass
    if len(nums) < 2:
        return None
    nums.sort()
    return {
        "n": len(nums),
        "min": nums[0],
        "q25": nums[len(nums) // 4],
        "median": statistics.median(nums),
        "q75": nums[(3 * len(nums)) // 4],
        "max": nums[-1],
    }


def load_rows(path: Path) -> list[dict]:
    import csv

    with path.open() as fh:
        return list(csv.DictReader(fh))


def resolve_manifest(results_dir: Path, target: str) -> Path:
    ranked = results_dir / f"ranked_{target}.csv"
    if ranked.is_file():
        return ranked
    struct = results_dir / f"structure_manifest_{target}.csv"
    if struct.is_file():
        return struct
    raise FileNotFoundError(f"no ranked_{target}.csv or structure_manifest_{target}.csv in {results_dir}")


def summarize(rows: list[dict], label: str, top_n: int) -> str:
    n = len(rows)
    keys = list(rows[0].keys()) if rows else []
    out = []
    p = out.append
    p(f"AMPBinderDesign panel summary: {label}  ({n} designs)")
    p("")

    cat_col = find_col(keys, "catalytic_ok", "catalytic_ok_bool")
    if cat_col:
        n_ok = sum(1 for r in rows if is_true(r.get(cat_col)))
        pct = 100.0 * n_ok / n if n else 0.0
        verdict = "STRONG" if pct >= 30 else "MARGINAL" if pct >= 10 else "WEAK"
        p(f"INHIBITOR GATE  catalytic_ok=True: {n_ok}/{n}  ({pct:.1f}%)   [{verdict}]")
    else:
        p("INHIBITOR GATE  catalytic_ok column not found (run active_site_overlap)")
    p("")

    iptm_col = find_col(keys, "boltz2_iptm")
    if cat_col and iptm_col:
        n_both = 0
        for r in rows:
            try:
                if is_true(r.get(cat_col)) and float(r.get(iptm_col)) >= 0.5:
                    n_both += 1
            except (TypeError, ValueError):
                pass
        p(f"STRUCTURE GATE  catalytic_ok and iPTM>=0.5: {n_both}/{n}  ({100.0 * n_both / n:.1f}%)")
        p("")

    gate_col = find_col(keys, "passes_structure_gates")
    if gate_col:
        n_gate = sum(1 for r in rows if is_true(r.get(gate_col)))
        p(f"STRUCTURE GATES  pass: {n_gate}/{n} ({100.0 * n_gate / n:.1f}%)")
        p("")

    p(f"{'metric':<42} {'n':>5} {'min':>9} {'q25':>9} {'median':>9} {'q75':>9} {'max':>9}")
    p("-" * 100)
    for raw, name in KEY_METRICS:
        col = find_col(keys, raw)
        if not col:
            continue
        s = numeric_stats([r.get(col, "") for r in rows])
        if s is None:
            continue
        p(
            f"{name:<42} {s['n']:>5} {s['min']:>9.3g} {s['q25']:>9.3g} "
            f"{s['median']:>9.3g} {s['q75']:>9.3g} {s['max']:>9.3g}"
        )
    p("")

    rank_col = find_col(keys, "rank_score")
    if rank_col:
        ranked = []
        for r in rows:
            try:
                ranked.append((float(r[rank_col]), r))
            except (TypeError, ValueError):
                pass
        ranked.sort(key=lambda x: -x[0])
        if ranked:
            cE = find_col(keys, "epitope_recall")
            cN = find_col(keys, "n_catalytic_contacts")
            cL = find_col(keys, "length")
            cQ = find_col(keys, "net_charge_pH7.4")
            cIP = find_col(keys, "boltz2_iptm", "iptm")
            cS = find_col(keys, "ipSAE_min")
            p(f"TOP {min(top_n, len(ranked))} BY rank_score:")
            p(
                f"  {'design_id':<28} {'rank':>7} {'iPTM':>6} {'ipSAE':>7} "
                f"{'eRec':>6} {'#cat':>5} {'L':>4} {'charge':>7}"
            )
            for sc, r in ranked[:top_n]:
                eid = (r.get("design_id", "") or "")[:28]
                er = r.get(cE, "") if cE else ""
                nc = r.get(cN, "") if cN else ""
                L = r.get(cL, "") if cL else ""
                ch = r.get(cQ, "") if cQ else ""
                ipt = r.get(cIP, "") if cIP else ""
                ipsae = r.get(cS, "") if cS else ""
                try:
                    er = f"{float(er):.2f}"
                except (TypeError, ValueError):
                    pass
                try:
                    ipt = f"{float(ipt):.2f}"
                except (TypeError, ValueError):
                    pass
                try:
                    ipsae = f"{float(ipsae):.3f}"
                except (TypeError, ValueError):
                    pass
                try:
                    ch = f"{float(ch):+.1f}"
                except (TypeError, ValueError):
                    pass
                p(f"  {eid:<28} {sc:>7.3f} {ipt:>6} {ipsae:>7} {er:>6} {nc:>5} {L:>4} {ch:>7}")
        p("")

    p("WET-LAB READINESS:")
    p("  Pool is already AMP-positive (AMPScanner + Macrel). Rank by Boltz-2 iPTM,")
    p("  ipSAE_min, and epitope coverage; require catalytic_ok for the shipped panel.")
    p("  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.")
    p("  - ipSAE_min is typically low here; treat it as a tie-breaker, not a hard cut.")
    p("  - Hemolysis is common among AMPs; prefer NonHemo but do not empty the panel.")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=None, help="single CSV (ranked or structure_manifest)")
    ap.add_argument("--results-dir", type=Path, default=REPO / "results")
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--out", default=None, help="write Markdown for a single --manifest")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="write summary_<TARGET>.md here when summarizing both targets",
    )
    args = ap.parse_args()

    if args.manifest:
        path = Path(args.manifest)
        rows = load_rows(path)
        if not rows:
            print("Empty manifest", file=sys.stderr)
            return 1
        text = summarize(rows, str(path), args.top_n)
        print(text)
        if args.out:
            outp = Path(args.out)
            outp.parent.mkdir(parents=True, exist_ok=True)
            outp.write_text(text + "\n")
            print(f"\nAlso wrote: {outp}")
        return 0

    out_dir = args.out_dir or (REPO / "results" / "plots" / "select")
    out_dir.mkdir(parents=True, exist_ok=True)
    for target in args.targets:
        path = resolve_manifest(args.results_dir, target)
        rows = load_rows(path)
        text = summarize(rows, f"{target}  {path}", args.top_n)
        print(text)
        print("")
        outp = out_dir / f"summary_{target}.md"
        outp.write_text(text + "\n")
        print(f"Wrote {outp}\n")
        sel = args.results_dir / f"selected_{target}.csv"
        if sel.is_file():
            sel_rows = load_rows(sel)
            if sel_rows:
                sel_text = summarize(sel_rows, f"{target} selected panel  {sel}", args.top_n)
                print(sel_text)
                print("")
                sel_out = out_dir / f"summary_selected_{target}.md"
                sel_out.write_text(sel_text + "\n")
                print(f"Wrote {sel_out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
