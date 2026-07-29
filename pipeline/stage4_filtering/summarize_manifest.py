#!/usr/bin/env python3
"""
Quality summary of a PepGen manifest CSV. Run after build_manifest / apply_filters
to decide whether the panel is wet-lab-ready.

Usage:
  python summarize_manifest.py --manifest results/manifest_<T>.csv [--top-n 10] [--out report.md]
"""
import argparse, csv, statistics, sys

KEY_METRICS = [
    ("epitope_recall",        "Active-site coverage (fraction)"),
    ("n_catalytic_contacts",  "# catalytic-core residues contacted"),
    ("interface_precision",   "Interface focus on active site"),
    ("complex_iptm",          "BoltzGen iPTM (binder-target interface)"),
    ("iptm",                  "iPTM"),
    ("length",                "Binder length (aa)"),
    ("net_charge_pH7.4",      "Net charge at pH 7.4"),
    ("hydrophobic_moment",    "Hydrophobic moment (amphipathicity)"),
    ("delivery_proxy",        "Periplasmic-delivery proxy [0..1]"),
    ("aggregation_proxy",     "Aggregation proxy (KD window max)"),
    ("n_liabilities",         "# synthesis liabilities"),
]


def find_col(keys, name):
    """Resolve a column name allowing build_manifest's 'tag.col' prefixing."""
    for k in keys:
        if k == name or k.endswith("." + name):
            return k
    return None


def numeric_stats(vals):
    nums = []
    for v in vals:
        try:
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


def fmt(x):
    if isinstance(x, float):
        return f"{x:.3g}"
    return str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--out", default=None, help="optional path to also write a Markdown report")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.manifest)))
    if not rows:
        print("Empty manifest", file=sys.stderr)
        return 1
    n = len(rows)
    keys = list(rows[0].keys())

    out = []
    p = out.append
    p(f"PepGen quality summary: {args.manifest}  ({n} designs)\n")

    # 1) the inhibitor gate: catalytic_ok yield
    cat_col = find_col(keys, "catalytic_ok")
    if cat_col:
        n_ok = sum(1 for r in rows if str(r.get(cat_col, "")).strip().lower() == "true")
        pct = 100.0 * n_ok / n
        verdict = ("STRONG"   if pct >= 30 else
                   "MARGINAL" if pct >= 10 else
                   "WEAK (epitope conditioning may have under-fired; consider tightening hotspots or re-running)")
        p(f"INHIBITOR GATE  catalytic_ok=True: {n_ok}/{n}  ({pct:.1f}%)   [{verdict}]")
    else:
        p("INHIBITOR GATE  catalytic_ok column not found (run Stage 3 active_site_overlap first)")
    p("")

    # 2) distributions over key metrics
    p(f"{'metric':<42} {'n':>5} {'min':>9} {'q25':>9} {'median':>9} {'q75':>9} {'max':>9}")
    p("-" * 100)
    for raw, label in KEY_METRICS:
        col = find_col(keys, raw)
        if not col:
            continue
        s = numeric_stats([r.get(col, "") for r in rows])
        if s is None:
            continue
        p(f"{label:<42} {s['n']:>5} {s['min']:>9.3g} {s['q25']:>9.3g} {s['median']:>9.3g} {s['q75']:>9.3g} {s['max']:>9.3g}")
    p("")

    # 3) gate pass rate
    pass_col = find_col(keys, "_pass")
    if pass_col:
        n_pass = sum(1 for r in rows if str(r.get(pass_col, "")).strip().lower() == "true")
        p(f"GATE PASS  {n_pass}/{n} ({100.0 * n_pass / n:.1f}%) survive all active filters")
        if n_pass < 32:
            p(f"  WARN: fewer than 32 designs pass; loosen filters.json (but never drop catalytic_ok).")
        p("")

    # 4) top N by rank_score
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
            cE  = find_col(keys, "epitope_recall")
            cN  = find_col(keys, "n_catalytic_contacts")
            cL  = find_col(keys, "length")
            cQ  = find_col(keys, "net_charge_pH7.4")
            cIP = find_col(keys, "complex_iptm") or find_col(keys, "iptm")
            p(f"TOP {min(args.top_n, len(ranked))} BY rank_score:")
            p(f"  {'design_id':<38} {'rank':>7} {'eRecall':>8} {'#cat':>5} {'iptm':>6} {'L':>4} {'charge':>7}")
            for sc, r in ranked[:args.top_n]:
                eid = (r.get("design_id", "") or "")[:38]
                er  = r.get(cE, "")  if cE  else ""
                nc  = r.get(cN, "")  if cN  else ""
                L   = r.get(cL, "")  if cL  else ""
                ch  = r.get(cQ, "")  if cQ  else ""
                ipt = r.get(cIP, "") if cIP else ""
                try: er  = f"{float(er):.2f}"
                except (TypeError, ValueError): pass
                try: ipt = f"{float(ipt):.2f}"
                except (TypeError, ValueError): pass
                try: ch  = f"{float(ch):+.1f}"
                except (TypeError, ValueError): pass
                p(f"  {eid:<38} {sc:>7.3f} {er:>8} {nc:>5} {ipt:>6} {L:>4} {ch:>7}")
        p("")

    # 5) brief decision guidance
    p("WET-LAB READINESS HEURISTIC:")
    p("  - catalytic_ok yield >= 30%:           strong; ship the top-32 panel + alanine controls.")
    p("  - catalytic_ok yield 10..30%:          marginal; ship but expect a lower hit rate in vitro.")
    p("  - catalytic_ok yield < 10%:            weak; rerun BoltzGen with a tighter hotspot list,")
    p("                                         or filter to catalytic_ok=True before selection.")
    p("  - median iPTM >= 0.7 in top-32:        good binder-target interface confidence.")
    p("  - median net charge in top-32 +2..+8:  reasonable periplasmic-delivery profile.")

    text = "\n".join(out)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write("```\n" + text + "\n```\n")
        print(f"\nAlso wrote: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
