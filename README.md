# Co2FeSn (Heusler) — QE inputs + lightweight analysis scripts

This repository is a **curated, public-friendly** snapshot of a DFT workflow for the Co₂FeSn Heusler alloy
using **Quantum ESPRESSO (pw.x)**, plus small Python utilities to summarize convergence series.

## Contents
- `structures/` — starting structure (`Co2FeSn_Prim.cif`)
- `inputs/` — example inputs + batch scripts (paths generalized)
  - `scf_Co2FeSn.in` — SCF input (edit `pseudo_dir` and pseudopotential filenames as needed)
  - `Optimizacion_ecut.sh`, `Optimizacion_kpoints.sh` — example scripts to scan `ecutwfc` / k-mesh
- `outputs_examples/` — **sanitized** example outputs (no personal/HPC paths)
- `scripts/` — parsing + summary scripts
- `results/` — CSV summaries of total energy vs `ecutwfc` and k-mesh

## Quick start
1. Install Quantum ESPRESSO and make sure `pw.x` is on your PATH.
2. Place the required pseudopotentials in `./pseudos/` (or change `pseudo_dir` in inputs).
3. Run an SCF:
   ```bash
   pw.x < inputs/scf_Co2FeSn.in > scf.out
   ```
4. Summarize the provided example series into CSV:
   ```bash
   python3 scripts/summarize_qe_series.py --root .
   ```

## Notes
- Any occurrences of cluster-specific paths were replaced with placeholders for portability.
- The `.out` files in `outputs_examples/` are meant as **examples**; for real runs you usually
  do not commit large scratch folders (`*.save/`, `outdir/`, etc.). See `.gitignore`.
