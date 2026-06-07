#!/usr/bin/env python3
"""
reproduce.py — reproduces every numerical claim in Version 4 of
"A Curvature Identity for the Fibonacci Digital-Root Cycle Modulo 9".

Pure standard library (math, cmath, itertools); no numpy required.
Run:  python3 reproduce.py
It prints each claim with its computed value and asserts the ones the paper
states as exact. A non-zero exit means a claim failed to reproduce.

Definitions (Section 2 of the paper):
  - Modulus m, Pisano period P = pi(m).
  - Orbit a_{n+1} = (a_n + a_{n-1}) mod m, restricted to seeds of minimal period P.
  - Display d(r) = m if r == 0 else r  (digital-root convention, values in 1..m).
  - Curvature kappa_n = d_{n+1} - 2 d_n + d_{n-1} (cyclic), S = sum kappa_n^2.
"""

import cmath
import math
import itertools


# ----- core machinery -------------------------------------------------------

def pisano(m):
    a, b = 0, 1
    for k in range(1, 6 * m + 1):
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            return k
    raise RuntimeError("Pisano period not found for m=%d" % m)


def disp(r, m):
    return m if r == 0 else r


def orbit_displayed(seed, m, P):
    a, b = seed
    out = []
    for _ in range(P):
        out.append(disp(a, m))
        a, b = b, (a + b) % m
    return out


def curvature_sum(d):
    N = len(d)
    return sum((d[(i + 1) % N] - 2 * d[i] + d[(i - 1) % N]) ** 2 for i in range(N))


def kappa_vector(d):
    N = len(d)
    return [d[(i + 1) % N] - 2 * d[i] + d[(i - 1) % N] for i in range(N)]


def full_period_orbits(m):
    """Yield one representative S per full-period orbit, with seed counts."""
    P = pisano(m)
    visited = set()
    orbits = []
    for a0 in range(m):
        for a1 in range(m):
            seed = (a0, a1)
            if seed in visited:
                continue
            pairs = []
            a, b = a0, a1
            for _ in range(P):
                pairs.append((a, b))
                a, b = b, (a + b) % m
            for p in pairs:
                visited.add(p)
            if len(set(pairs)) != P or (a, b) != seed:
                continue  # not minimal period P
            d = [disp(p[0], m) for p in pairs]
            orbits.append({"seed": seed, "S": curvature_sum(d), "pairs": pairs})
    return P, orbits


def smax(m):
    _, orbits = full_period_orbits(m)
    return max(o["S"] for o in orbits)


def fib_S(m):
    P = pisano(m)
    return curvature_sum(orbit_displayed((1, 1), m, P))


# ----- spectral identity (Proposition 1) ------------------------------------

def spectral_S(d):
    N = len(d)
    X = [sum(d[n] * cmath.exp(-2j * math.pi * k * n / N) for n in range(N)) for k in range(N)]
    mu2 = [16 * math.sin(math.pi * k / N) ** 4 for k in range(N)]
    return sum(mu2[k] * abs(X[k]) ** 2 for k in range(N)) / N


def near_nyquist_fraction(d):
    """Fraction of S in modes with sin^2(pi k/N) >= 1/2 (the near-Nyquist half)."""
    N = len(d)
    X = [sum(d[n] * cmath.exp(-2j * math.pi * k * n / N) for n in range(N)) for k in range(N)]
    c = [16 * math.sin(math.pi * k / N) ** 4 * abs(X[k]) ** 2 / N for k in range(N)]
    S = sum(c)
    band = sum(c[k] for k in range(N) if math.sin(math.pi * k / N) ** 2 >= 0.5 - 1e-12)
    return 100 * band / S


# ----- checks ---------------------------------------------------------------

def main():
    print("== Theorem 1: trichotomy modulo 9 ==")
    P, orbits = full_period_orbits(9)
    from collections import Counter
    by_seed = Counter()
    for o in orbits:
        by_seed[o["S"]] += len(o["pairs"])
    n_full = sum(by_seed.values())
    print("  P(9) =", P, "| full-period seeds =", n_full,
          "| distribution =", dict(sorted(by_seed.items())))
    assert P == 24 and n_full == 72
    assert dict(by_seed) == {648: 24, 972: 24, 1296: 24}
    print("  S((1,1)) =", fib_S(9), "= max")
    assert fib_S(9) == 1296 == smax(9)

    print("== kappa vector and half sums (Table 1, Section 4.1) ==")
    d = orbit_displayed((1, 1), 9, 24)
    kap = kappa_vector(d)
    print("  kappa =", kap)
    print("  sum kappa =", sum(kap), "| sum kappa^2 =", sum(k * k for k in kap))
    h1 = sum(k * k for k in kap[:12])
    h2 = sum(k * k for k in kap[12:])
    print("  half sums:", h1, "and", h2, "(imbalance %d)" % (h2 - h1))
    assert sum(kap) == 0 and sum(k * k for k in kap) == 1296
    assert (h1, h2) == (459, 837)
    breaks = [n for n in range(12) if kap[n + 12] != -kap[n]]
    print("  antisymmetry kappa_{n+12} = -kappa_n breaks at n =", breaks)
    assert breaks == [0, 10, 11]

    print("== Remark: three values are the three unit classes (Z/9)*/{+-1} ==")
    res9 = [1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 0, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1, 0]
    scaled = {g: curvature_sum([disp((g * r) % 9, 9) for r in res9]) for g in (1, 2, 4, 5, 7, 8)}
    print("  S(g * orbit):", scaled)
    assert scaled[1] == scaled[8] == 1296
    assert scaled[2] == scaled[7] == 972
    assert scaled[4] == scaled[5] == 648

    print("== Proposition 1: spectral identity reproduces S exactly ==")
    for m, seed in [(9, (1, 1)), (9, (1, 2)), (11, (1, 1)), (16, (1, 1)), (27, (1, 1))]:
        dd = orbit_displayed(seed, m, pisano(m))
        direct, spec = curvature_sum(dd), spectral_S(dd)
        print("  m=%2d seed=%s  S_direct=%d  S_spectral=%.6f" % (m, seed, direct, spec.real))
        assert abs(direct - spec.real) < 1e-6

    print("== near-Nyquist half fractions (Section 4, Remark, Section 5) ==")
    for g in (1, 2, 4):
        dd = [disp((g * r) % 9, 9) for r in res9]
        print("  m=9 class S=%d : %.1f%%" % (curvature_sum(dd), near_nyquist_fraction(dd)))
    for seed in [(0, 3), (1, 1)]:
        dd = orbit_displayed(seed, 40, pisano(40))
        print("  m=40 seed=%s S=%d : %.1f%%" % (seed, curvature_sum(dd), near_nyquist_fraction(dd)))

    print("== Section 5: maximality sweep, m <= 120 with even Pisano period ==")
    wins, fails = [], []
    for m in range(2, 121):
        if pisano(m) % 2:
            continue
        (wins if fib_S(m) == smax(m) else fails).append(m)
    print("  even-Pisano moduli =", len(wins) + len(fails),
          "| (1,1) wins =", len(wins), "| fails =", len(fails))
    print("  first failures:", fails[:3], "| full failure list:", fails)
    assert len(wins) + len(fails) == 118 and len(wins) == 78 and len(fails) == 40
    assert fails[:3] == [15, 33, 40]

    print("== Section 6: S_max = m*(pi/2)^2 is a two-point coincidence ==")
    for m in [3, 9, 27, 81, 243]:
        P = pisano(m)
        S = smax(m)
        pred = m * (P // 2) ** 2
        print("  m=%3d  S_max=%d  pred=%d  %s" % (m, S, pred, "HOLDS" if S == pred else "fails"))
        if m in (9, 27):
            assert S == pred
        else:
            assert S != pred

    print("\nAll Version 4 claims reproduced. OK.")


if __name__ == "__main__":
    main()
