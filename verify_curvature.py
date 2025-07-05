#!/usr/bin/env python3
"""
verify_curvature.py

Exhaustively verifies the curvature trichotomy {648, 972, 1296}
for all 81 ordered seeds (a0, a1) with digital-root values 1-9.
"""

import itertools
from collections import Counter

def to_residue(x: int) -> int:
    """Map 9 → 0 so that arithmetic is clean mod 9."""
    return 0 if x == 9 else x

def display_val(r: int) -> int:
    """Map 0 → 9 for digital-root display."""
    return 9 if r == 0 else r

def orbit(seed):
    """Return (displayed orbit, raw residues) of length 24."""
    r_prev, r_curr = map(to_residue, seed)
    disp = [display_val(r_prev), display_val(r_curr)]
    resid = [r_prev, r_curr]
    for _ in range(22):           # need 24 terms in total
        r_next = (r_prev + r_curr) % 9
        resid.append(r_next)
        disp.append(display_val(r_next))
        r_prev, r_curr = r_curr, r_next
    return disp, resid

def minimal_period(res):
    """Smallest n ≥ 1 with (a_n, a_{n+1}) = (a_0, a_1)."""
    r0, r1 = res[0], res[1]
    r_prev, r_curr = r0, r1
    for n in range(1, 25):
        r_next = (r_prev + r_curr) % 9
        r_prev, r_curr = r_curr, r_next
        if (r_prev, r_curr) == (r0, r1):
            return n
    return 24

def curvature_sum(disp):
    """Σ κ_n² with κ_n = a_{n+1} − 2a_n + a_{n−1} (indices mod 24)."""
    return sum((disp[(i+1) % 24] - 2*disp[i] + disp[(i-1) % 24])**2
               for i in range(24))

def main():
    counts, full = Counter(), []
    for a0, a1 in itertools.product(range(1, 10), repeat=2):
        disp, res = orbit((a0, a1))
        if minimal_period(res) == 24:
            S = curvature_sum(disp)
            counts[S] += 1
            full.append(((a0, a1), S))

    print("Full-period seeds:", sum(counts.values()))   # should be 72
    print("Curvature distribution:", dict(counts))      # {648:24, 972:24, 1296:24}
    print("Seeds with S = 1296:",
          [seed for seed, S in full if S == 1296])

    assert counts == Counter({648: 24, 972: 24, 1296: 24}), \
        "Distribution does not match theorem."

if __name__ == "__main__":
    main()