#!/usr/bin/env python3
"""Small helper functions to parse Quantum ESPRESSO (pw.x) output files.

This is intentionally lightweight and dependency-free (standard library only).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

_RE_TOTAL_ENERGY = re.compile(r"!\s+total energy\s+=\s+([\-\d\.Ee\+]+)\s+Ry")

def parse_total_energy_ry(text: str) -> Optional[float]:
    """Return the *last* reported total energy (Ry) or None if not found."""
    matches = _RE_TOTAL_ENERGY.findall(text)
    if matches:
        return float(matches[-1])
    # fallback: sometimes 'total energy' appears without '!'
    matches = re.findall(r"total energy\s+=\s+([\-\d\.Ee\+]+)\s+Ry", text)
    if matches:
        return float(matches[-1])
    return None

def parse_total_energy_from_file(path: Path) -> Optional[float]:
    return parse_total_energy_ry(path.read_text(errors="ignore"))

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Parse last total energy (Ry) from QE pw.x output.")
    ap.add_argument("file", type=Path, help="pw.x output file (e.g., scf.out)")
    args = ap.parse_args()
    e = parse_total_energy_from_file(args.file)
    if e is None:
        raise SystemExit("No total energy found.")
    print(e)
