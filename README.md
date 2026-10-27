# A Curvature Identity for the Fibonacci Digital-Root Cycle Modulo 9

This repository contains the source files and compiled PDFs for the paper:

**A Curvature Identity for the Fibonacci Digital-Root Cycle Modulo 9**
*Guy Rinzema*

**DOI (all versions):** [10.5281/zenodo.15811767](https://doi.org/10.5281/zenodo.15811767)

## Overview

For each seed we form the 24-periodic orbit of the Fibonacci recurrence modulo 9, write
its terms with the digital-root convention (0 displayed as 9), and measure its *digital
curvature* — the sum of squared discrete second differences. The paper makes two points:

1. **Trichotomy (modulo 9).** Over the 72 full-period seeds the curvature takes exactly
   three values — 648, 972, 1296 — each on 24 seeds. The maximum 1296 belongs precisely to
   the classical Fibonacci cycle (1,1) and its rotations. This is proved by exhaustive
   enumeration. Modulo 9 the 0 → 9 display rule does not matter: raw residues 0–8 give the
   same three values, because every orbit negates itself after half a period.
2. **Spectral reformulation.** Curvature is the squared norm of the cyclic second-difference
   (circulant) operator, so it has an exact Fourier form, S = (16/N) Σ sin⁴(πk/N) |d̂ₖ|².
   It is a high-pass energy, and modulo 9 the Fibonacci class is the one with the most of
   its curvature near Nyquist (≈98%). That describes the maximizer but does not explain it:
   at other moduli the near-Nyquist energy does not identify the maximizing orbit. The
   maximality itself is a small-modulus phenomenon (it fails for 40 of the 118 moduli
   ≤ 120 with even Pisano period, first at m = 15, 33, 40), and the earlier conjecture
   S_max = m(π(m)/2)² holds only at m = 9 and 27 among the moduli ≤ 300 checked.

## Repository layout

- `papers/v5/` — **current version** (source and PDF). Corrects errors in V4 found by an
  independent fact-check; see the revision history in the paper. Includes `reproduce.py`,
  which asserts every numerical claim in V5 over the stated ranges. On Zenodo as
  [10.5281/zenodo.23687697](https://doi.org/10.5281/zenodo.23687697) (Version 5 DOI).
  **Supersedes V1–V4.**
- `papers/v4/` — Version 4 (source and PDF). Superseded: its claims that the curvature
  maximizer is always the most near-Nyquist orbit, that raw residues give a different S
  modulo 9, and that every claim was verified by two independent implementations are
  false or overstated.
- `papers/v3/`, `papers/v2/`, `papers/v1/` — earlier versions, retained as the historical
  record. **Note:** these contain a non-rigorous "proof" of the trichotomy, a false
  half-period symmetry argument (V2), and the false conjecture S_max = m(π(m)/2)²; all are
  corrected in V4 and V5.
- `verify_curvature.py` — the original modulo-9 trichotomy check that accompanies V1–V3.

## Reproducing the results

```
python3 papers/v5/reproduce.py                # asserts every Version 5 claim
python3 verify_curvature.py                   # original modulo-9 trichotomy check
```

All use only the Python standard library and run in seconds.

## License

This work is made available under the **Creative Commons Attribution 4.0 International
License (CC BY 4.0)**.

## Author

Guy Rinzema
[https://github.com/founder-guy](https://github.com/founder-guy)
