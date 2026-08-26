#!/usr/bin/env python3
"""Batch ipSAE for Boltz-2 complexes using Dunbrack's ipsae.py.

For each prediction folder, runs:
  python ipsae.py pae_<stem>_model_0.npz <stem>_model_0.cif <pae_cutoff> <dist_cutoff>

ipSAE is asymmetric (A→B ≠ B→A). This reports ipSAE_min of the two Type=asym
rows for chains A/B, matching the IDPepDesign / Forge convention.

Usage (from AMPBinderDesign, after boltz predict --write_full_pae):
  python src/4_structure_prediction/run_ipsae.py
  python src/4_structure_prediction/run_ipsae.py --boltz-out boltz_results --workers 8

Writes results/ipsae_NDM5.csv and results/ipsae_KPC3.csv for build_manifest.py.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
IPSAE_SCRIPT = HERE / "ipsae.py"
TARGETS = ("NDM5", "KPC3")
TARGET_CHAIN = "A"
BINDER_CHAIN = "B"
IPSAE_EXTRA = ("ipSAE_d0chn", "ipSAE_d0dom", "ipTM_af", "pDockQ", "pDockQ2", "LIS")
_RANK_PREFIX_RE = re.compile(r"^rank\d+_", re.IGNORECASE)
_MODEL_SUFFIX_RE = re.compile(r"_model_\d+$", re.IGNORECASE)
_SKIP_DIR_PARTS = {"processed", "mols", "msa"}

COLS = [
    "design_id",
    "ipSAE_min",
    "ipSAE_max",
    "ipSAE_d0chn",
    "ipSAE_d0dom",
    "ipTM_af",
    "pDockQ",
    "pDockQ2",
    "LIS",
    "status",
]


def design_id_from_stem(stem: str) -> str:
    stem = _RANK_PREFIX_RE.sub("", stem)
    return _MODEL_SUFFIX_RE.sub("", stem)


def as_float(row: dict, key: str):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return None


def parse_ipsae_summary(path: Path, target_chain: str, binder_chain: str) -> dict:
    if not path.exists():
        return {}
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    header_idx = next((i for i, line in enumerate(lines) if line.startswith("Chn1,")), None)
    if header_idx is None:
        return {}
    header = [p.strip() for p in lines[header_idx].split(",")]
    rows = []
    for line in lines[header_idx + 1 :]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))

    wanted = {target_chain, binder_chain}
    rows = [r for r in rows if {r.get("Chn1"), r.get("Chn2")} == wanted] or rows
    if not rows:
        return {}

    asym = [v for v in (as_float(r, "ipSAE") for r in rows if r.get("Type") == "asym") if v is not None]
    max_rows = [r for r in rows if r.get("Type") == "max"]
    base = max_rows[0] if max_rows else max(rows, key=lambda r: as_float(r, "ipSAE") or 0.0)

    out = {}
    if asym:
        out["ipSAE_min"] = min(asym)
        out["ipSAE_max"] = max(asym)
    else:
        v = as_float(base, "ipSAE")
        if v is not None:
            out["ipSAE_max"] = v
    for key in IPSAE_EXTRA:
        v = as_float(base, key)
        if v is not None:
            out[key] = v
    return out


def ipsae_txt_path(cif: Path, pae_cutoff: float, dist_cutoff: float) -> Path:
    pae_s = f"{int(pae_cutoff):02d}" if pae_cutoff < 10 else str(int(pae_cutoff))
    dist_s = f"{int(dist_cutoff):02d}" if dist_cutoff < 10 else str(int(dist_cutoff))
    return cif.with_name(f"{cif.stem}_{pae_s}_{dist_s}.txt")


def run_ipsae(pae: Path, cif: Path, pae_cutoff: float, dist_cutoff: float) -> Path:
    out_txt = ipsae_txt_path(cif, pae_cutoff, dist_cutoff)
    if out_txt.exists() and out_txt.stat().st_size > 0:
        return out_txt
    cmd = [
        sys.executable,
        str(IPSAE_SCRIPT),
        str(pae.resolve()),
        str(cif.resolve()),
        str(int(pae_cutoff) if pae_cutoff == int(pae_cutoff) else pae_cutoff),
        str(int(dist_cutoff) if dist_cutoff == int(dist_cutoff) else dist_cutoff),
    ]
    proc = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(err[-800:] or f"ipsae.py exited {proc.returncode}")
    if not out_txt.exists():
        raise FileNotFoundError(f"ipsae.py did not write {out_txt}")
    return out_txt


def discover_jobs(boltz_out: Path, target: str | None) -> list[tuple[str, Path | None, Path]]:
    jobs = []
    for cif in sorted(boltz_out.rglob("*_model_0.cif")):
        if any(part in _SKIP_DIR_PARTS for part in cif.parts):
            continue
        did = design_id_from_stem(cif.stem)
        if target and target not in did:
            continue
        pae = cif.with_name(f"pae_{cif.stem}.npz")
        if not pae.exists():
            hits = sorted(cif.parent.glob("pae_*.npz"))
            pae = hits[0] if hits else None
        jobs.append((did, pae, cif))
    return jobs


def score_one(did: str, pae_s: str | None, cif_s: str, pae_cutoff: float, dist_cutoff: float) -> dict:
    row = {c: "" for c in COLS}
    row["design_id"] = did
    if not pae_s:
        row["status"] = "missing_pae"
        return row
    pae, cif = Path(pae_s), Path(cif_s)
    try:
        out_txt = run_ipsae(pae, cif, pae_cutoff, dist_cutoff)
        parsed = parse_ipsae_summary(out_txt, TARGET_CHAIN, BINDER_CHAIN)
        if not parsed:
            row["status"] = "ipsae_parse_failed"
            return row
        for k, v in parsed.items():
            if k in row:
                row[k] = v
        row["status"] = "ok"
    except Exception as exc:
        row["status"] = f"ERROR {type(exc).__name__}: {exc}"
    return row


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLS})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--boltz-out",
        default=os.path.join(REPO, "boltz_results"),
        help="Boltz-2 output dir (searched recursively for *_model_0.cif + pae_*.npz)",
    )
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--pae-cutoff", type=float, default=10)
    ap.add_argument("--dist-cutoff", type=float, default=10)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument(
        "--out-dir",
        default=os.path.join(REPO, "results"),
        help="writes ipsae_<TARGET>.csv here",
    )
    args = ap.parse_args()

    if not IPSAE_SCRIPT.is_file():
        raise FileNotFoundError(IPSAE_SCRIPT)

    boltz_out = Path(args.boltz_out)
    if not boltz_out.exists():
        print(f"Boltz output not found: {boltz_out}", file=sys.stderr)
        return 1

    for target in args.targets:
        jobs = discover_jobs(boltz_out, target)
        if not jobs:
            print(f"[{target}] no *_model_0.cif under {boltz_out}")
            continue
        rows: list[dict] = []
        n_ok = n_miss = n_fail = 0
        payloads = [
            (did, str(pae) if pae else None, str(cif), args.pae_cutoff, args.dist_cutoff)
            for did, pae, cif in jobs
        ]
        if args.workers > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futs = [pool.submit(score_one, *p) for p in payloads]
                for fut in as_completed(futs):
                    row = fut.result()
                    rows.append(row)
                    st = row.get("status", "")
                    if st == "ok":
                        n_ok += 1
                    elif st == "missing_pae":
                        n_miss += 1
                    else:
                        n_fail += 1
            rows.sort(key=lambda r: r["design_id"])
        else:
            for p in payloads:
                row = score_one(*p)
                rows.append(row)
                st = row.get("status", "")
                if st == "ok":
                    n_ok += 1
                elif st == "missing_pae":
                    n_miss += 1
                else:
                    n_fail += 1
        out = Path(args.out_dir) / f"ipsae_{target}.csv"
        write_csv(rows, out)
        print(
            f"[{target}] ipSAE ok={n_ok} missing_pae={n_miss} failed={n_fail} "
            f"/ {len(rows)} -> {os.path.relpath(out, REPO)}"
        )
        if n_miss and n_ok == 0:
            print(
                "  No PAE npz files. Re-run boltz with --write_full_pae "
                "(expected pae_<stem>_model_0.npz next to the cif).",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
