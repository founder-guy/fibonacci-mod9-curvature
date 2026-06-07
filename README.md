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
   enumeration.
2. **Spectral reformulation.** Curvature is the squared norm of the cyclic second-difference
   (circulant) operator, so it has an exact Fourier form, S = (16/N) Σ sin⁴(πk/N) |d̂ₖ|².
   It is a high-pass energy: the Fibonacci orbit maximizes it because it is the most
   high-frequency of the three classes (≈98% of its energy is near Nyquist). The same lens
   shows the maximality is a small-modulus phenomenon (it fails for 40 of the 118 moduli
   ≤ 120 with even Pisano period, first at m = 15, 33, 40) and retires the earlier
   conjecture S_max = m(π(m)/2)², which holds only at m = 9 and 27.

## Versions

- `papers/v4/` — **current version.** Spectral reformulation; corrects errors in V1–V3
  (see the revision history / errata in the paper). Includes `reproduce.py`, which
  reproduces and asserts every numerical claim. **Supersedes V1–V3.**
- `papers/v3/`, `papers/v2/`, `papers/v1/` — earlier versions, retained as the historical
  record. **Note:** these contain a non-rigorous "proof" of the trichotomy, a false
  half-period symmetry argument (V2), and the false conjecture S_max = m(π(m)/2)²; all are
  corrected in V4.

## Reproducing the results

```
python3 papers/v4/reproduce.py    # reproduces and asserts every Version 4 claim
python3 verify_curvature.py       # original modulo-9 trichotomy check
```

Both use only the Python standard library.

## License

This work is made available under the **Creative Commons Attribution 4.0 International
License (CC BY 4.0)**.

## Author

Guy Rinzema
[https://github.com/founder-guy](https://github.com/founder-guy)
