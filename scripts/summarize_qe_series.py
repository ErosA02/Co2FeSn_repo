#!/usr/bin/env python3
"""Summarize convergence series (ecut / k-mesh) from QE outputs.

Expected layout (can be changed via CLI args):
- outputs_examples/opt_ecut/*.out  -> results/energy_vs_ecut.csv
- outputs_examples/opt_kpoints/*.out -> results/energy_vs_kmesh.csv
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import List, Dict

from parse_qe_output import parse_total_energy_from_file

def summarize_ecut(files: List[Path]) -> List[Dict[str, object]]:
    rows = []
    for f in files:
        m = re.search(r"_(\d+)\.out$", f.name)
        if not m:
            continue
        ecut = int(m.group(1))
        e = parse_total_energy_from_file(f)
        rows.append({"ecutwfc_Ry": ecut, "total_energy_Ry": e, "file": str(f)})
    return sorted(rows, key=lambda r: r["ecutwfc_Ry"])

def summarize_kmesh(files: List[Path]) -> List[Dict[str, object]]:
    rows = []
    for f in files:
        m = re.search(r"_(\d+)_kpoints\.out$", f.name)
        if not m:
            continue
        k = int(m.group(1))
        e = parse_total_energy_from_file(f)
        rows.append({"k_mesh": k, "total_energy_Ry": e, "file": str(f)})
    return sorted(rows, key=lambda r: r["k_mesh"])

def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."), help="Repository root (default: .)")
    args = ap.parse_args()

    ecut_dir = args.root / "outputs_examples" / "opt_ecut"
    k_dir = args.root / "outputs_examples" / "opt_kpoints"

    ecut_files = sorted(ecut_dir.glob("*.out"))
    k_files = sorted(k_dir.glob("*.out"))

    ecut_rows = summarize_ecut(ecut_files)
    k_rows = summarize_kmesh(k_files)

    write_csv(args.root / "results" / "energy_vs_ecut.csv", ecut_rows)
    write_csv(args.root / "results" / "energy_vs_kmesh.csv", k_rows)

    print(f"Wrote {len(ecut_rows)} rows -> results/energy_vs_ecut.csv")
    print(f"Wrote {len(k_rows)} rows -> results/energy_vs_kmesh.csv")

if __name__ == "__main__":
    main()
