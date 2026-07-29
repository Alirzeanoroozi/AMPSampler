#!/usr/bin/env python3
"""
Stage 3 - Parse Boltz-2 predictions into a design_id-keyed score table.

Boltz-2 writes, per prediction <name>:
  predictions/<name>/confidence_<name>_model_0.json   (iptm, ptm, complex_plddt, ...)
  predictions/<name>/affinity_<name>.json             (affinity_pred_value, affinity_probability_binary)
The Boltz-2 affinity head is used as ONE ranking signal (not the only one). `<name>` must
be the design_id so scores join to the manifest.

Output CSV columns: design_id, boltz2_iptm, boltz2_ptm, boltz2_plddt,
                     boltz2_affinity_pred_value, boltz2_affinity_prob_binary
Runs in the base conda env. Usage:
  python parse_boltz2.py --boltz_out <dir> --out boltz2_<T>.csv
"""
import argparse, csv, glob, json, os


def find(d, *keys):
    for k in keys:
        if k in d:
            return d[k]
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz_out", required=True, help="Boltz-2 output dir (contains predictions/)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    conf_files = glob.glob(os.path.join(args.boltz_out, "**", "confidence_*.json"), recursive=True)
    rows = []
    for cf in sorted(conf_files):
        name = os.path.basename(os.path.dirname(cf))
        try:
            conf = json.load(open(cf))
        except Exception:
            continue
        aff_path = glob.glob(os.path.join(os.path.dirname(cf), "affinity_*.json"))
        aff = {}
        if aff_path:
            try:
                aff = json.load(open(aff_path[0]))
            except Exception:
                aff = {}
        rows.append({
            "design_id": name,
            "boltz2_iptm": find(conf, "iptm", "complex_iptm"),
            "boltz2_ptm": find(conf, "ptm"),
            "boltz2_plddt": find(conf, "complex_plddt", "plddt"),
            "boltz2_affinity_pred_value": find(aff, "affinity_pred_value"),
            "boltz2_affinity_prob_binary": find(aff, "affinity_probability_binary"),
        })

    cols = ["design_id", "boltz2_iptm", "boltz2_ptm", "boltz2_plddt",
            "boltz2_affinity_pred_value", "boltz2_affinity_prob_binary"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"Parsed {len(rows)} Boltz-2 predictions -> {args.out}")


if __name__ == "__main__":
    main()
