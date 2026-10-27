#!/usr/bin/env python3
"""
reproduce.py -- checks the numerical claims in Version 5 of
"A Curvature Identity for the Fibonacci Digital-Root Cycle Modulo 9".

Standard library only. Run:  python3 papers/v5/reproduce.py
Exact claims are asserted exactly. Finite-range claims are asserted over the
range the paper states. Spectral quantities are floating point; every asserted
comparison has a margin far larger than rounding error. A non-zero exit means a
claim failed to reproduce.

Definitions (Section 2 of the paper):
  - Modulus m >= 2, Pisano period P = pi(m).
  - Orbit a_{n+1} = a_n + a_{n-1} (mod m); "full-period" means minimal period P.
  - Raw residues r in {0..m-1}; display d(r) = m if r == 0 else r.
  - kappa_n = d_{n+1} - 2 d_n + d_{n-1} (cyclic), S = sum kappa_n^2.
  - Near-Nyquist half: modes k with P/4 <= k <= 3P/4. B = curvature energy there.
"""

import math
from collections import Counter
from math import gcd


# ----- orbits ---------------------------------------------------------------

def orbits(m, p=1, q=1):
    """All cycles of (a, b) -> (b, p*b + q*a mod m), each a list of pairs.

    The map must be a bijection: always for p = q = 1, and for q = 2 when m is odd.
    """
    seen, out = set(), []
    for a in range(m):
        for b in range(m):
            if (a, b) in seen:
                continue
            cyc, x = [], (a, b)
            while x not in seen:
                seen.add(x)
                cyc.append(x)
                x = (x[1], (p * x[1] + q * x[0]) % m)
            assert x == (a, b), "map is not a bijection mod %d" % m
            out.append(cyc)
    return out


def orbit_of(cycles, pair):
    return next(c for c in cycles if pair in c)


def starting_at(cyc, pair):
    i = cyc.index(pair)
    return cyc[i:] + cyc[:i]


def shown(cyc, m):
    return [a if a else m for a, _ in cyc]


def raw(cyc):
    return [a for a, _ in cyc]


def kappa(x):
    N = len(x)
    return [x[(n + 1) % N] - 2 * x[n] + x[n - 1] for n in range(N)]


def S(x):
    return sum(v * v for v in kappa(x))


def matpow_is_minus_identity(m, e):
    """Is M^e == -I (mod m) for M = [[0,1],[1,1]]?"""
    R, A = [[1, 0], [0, 1]], [[0, 1], [1, 1]]
    while e:
        if e & 1:
            R = [[sum(R[i][t] * A[t][j] for t in range(2)) % m for j in range(2)] for i in range(2)]
        A = [[sum(A[i][t] * A[t][j] for t in range(2)) % m for j in range(2)] for i in range(2)]
        e >>= 1
    return R == [[m - 1, 0], [0, m - 1]]


# ----- spectra --------------------------------------------------------------

_twiddles = {}


def twiddles(N):
    if N not in _twiddles:
        _twiddles[N] = [complex(math.cos(2 * math.pi * j / N), -math.sin(2 * math.pi * j / N))
                        for j in range(N)]
    return _twiddles[N]


def dft(x, k):
    w, N = twiddles(len(x)), len(x)
    return sum(x[n] * w[(k * n) % N] for n in range(N))


def spectral_S(x):
    """Proposition 2 as stated: (16/N) sum_k sin^4(pi k/N) |x_hat_k|^2."""
    N = len(x)
    return sum(16 * math.sin(math.pi * k / N) ** 4 * abs(dft(x, k)) ** 2 for k in range(N)) / N


def mode_energy(x, k):
    """c_k = |DFT(kappa)_k|^2 / N, which equals (16/N) sin^4(pi k/N) |x_hat_k|^2."""
    return abs(dft(kappa(x), k)) ** 2 / len(x)


def near_nyquist(N):
    return [k for k in range(N) if N <= 4 * k <= 3 * N]


def B(x):
    K, N = kappa(x), len(x)
    return sum(abs(dft(K, k)) ** 2 for k in near_nyquist(N)) / N


# ----- checks ---------------------------------------------------------------

def main():
    print("== Theorem 1: trichotomy modulo 9, seed-period split ==")
    cy9 = orbits(9)
    assert sorted(len(c) for c in cy9) == [1, 8, 24, 24, 24]      # 81 = 1 + 8 + 72
    eight = next(c for c in cy9 if len(c) == 8)
    assert all(a % 3 == 0 and b % 3 == 0 for a, b in eight)        # period-8 seeds: both entries divisible by 3
    by_value = Counter()
    for c in cy9:
        if len(c) == 24:
            by_value[S(shown(c, 9))] += len(c)
    print("  seeds by period:", dict(Counter(len(c) for c in cy9 for _ in c)),
          "| curvature on full-period seeds:", dict(sorted(by_value.items())))
    assert dict(by_value) == {648: 24, 972: 24, 1296: 24}
    p9 = starting_at(orbit_of(cy9, (1, 1)), (1, 1))
    d9, r9 = shown(p9, 9), raw(p9)
    assert d9 == [1, 1, 2, 3, 5, 8, 4, 3, 7, 1, 8, 9, 8, 8, 7, 6, 4, 1, 5, 6, 2, 8, 1, 9]
    assert S(d9) == 1296 == max(S(shown(c, 9)) for c in cy9 if len(c) == 24)

    print("== Remark 1: representatives ==")
    p11 = orbit_of(orbits(11), (1, 1))
    sym9 = [a if a <= 4 else a - 9 for a in r9]
    print("  m=11 (1,1): raw S=%d, displayed S=%d | m=9 symmetric reps S=%d"
          % (S(raw(p11)), S(shown(p11, 11)), S(sym9)))
    assert (S(raw(p11)), S(shown(p11, 11))) == (638, 1430)
    assert S(sym9) == 648

    print("== Proposition 1: half-period negation ==")
    assert matpow_is_minus_identity(9, 12)
    for c in cy9:
        assert S(raw(c)) == S(shown(c, 9))                          # every orbit mod 9
        if len(c) == 24:
            rc, dc = raw(c), shown(c, 9)
            assert all(dc[n] == 9 - rc[(n + 12) % 24] for n in range(24))
            assert all(abs(mode_energy(dc, k) - mode_energy(rc, k)) < 1e-6 for k in range(24))
    halfneg, equal = [], []
    for m in range(2, 301):
        pr = orbit_of(orbits(m), (0, 1))
        P = len(pr)
        if P % 2:
            continue
        if matpow_is_minus_identity(m, P // 2):
            halfneg.append(m)
        if S(raw(pr)) == S(shown(pr, m)):
            equal.append(m)
    upto120 = [m for m in halfneg if m <= 120]
    print("  M^{P/2} == -I for %d of the even-period moduli <= 120:" % len(upto120), upto120)
    assert len(upto120) == 43
    assert halfneg == equal                                          # no mismatch for m <= 300
    print("  equal raw/displayed curvature <=> M^{P/2} == -I for every even-period m <= 300")

    print("== Remark 3: unit classes ==")
    scaled = {g: S([((g * a) % 9) or 9 for a in r9]) for g in (1, 2, 4, 5, 7, 8)}
    print("  S(g * orbit):", scaled)
    assert scaled[1] == scaled[8] == 1296
    assert scaled[2] == scaled[7] == 972
    assert scaled[4] == scaled[5] == 648

    print("== Table 1 and Section 4.1: second differences and half sums ==")
    kd, kr, ks = kappa(d9), kappa(r9), kappa(sym9)
    halves = lambda k: (sum(v * v for v in k[:12]), sum(v * v for v in k[12:]))
    breaks = lambda k: [n for n in range(12) if k[n + 12] != -k[n]]
    assert kd == [8, 1, 0, 1, 1, -7, 3, 5, -10, 13, -6, -2,
                  1, -1, 0, -1, -1, 7, -3, -5, 10, -13, 15, -16]
    assert sum(kd) == 0
    print("  displayed halves %s breaks %s | raw halves %s breaks %s | symmetric halves %s breaks %s"
          % (halves(kd), breaks(kd), halves(kr), breaks(kr), halves(ks), breaks(ks)))
    assert halves(kd) == (459, 837) and breaks(kd) == [0, 10, 11]
    assert halves(kr) == (837, 459) and breaks(kr) == [0, 10, 11]
    assert halves(ks) == (324, 324) and breaks(ks) == []

    print("== Proposition 2: spectral identity ==")
    for m, seed in [(9, (1, 1)), (9, (1, 2)), (11, (1, 1)), (16, (1, 1)), (27, (1, 1))]:
        x = shown(orbit_of(orbits(m), seed), m)
        print("  m=%2d seed=%s  S=%d  spectral=%.6f" % (m, seed, S(x), spectral_S(x)))
        assert abs(S(x) - spectral_S(x)) < 1e-6

    print("== Table 2 and Remark 4: spectrum of the (1,1) orbit, near-Nyquist shares ==")
    c = [mode_energy(d9, k) for k in range(13)]
    power = [abs(dft(d9, k)) ** 2 / 24 for k in range(13)]
    table_c = [0.000, 0.083, 0.242, 6.563, 3.375, 5.122, 13.500, 14.774, 30.375,
               222.937, 47.008, 277.022, 54.000]
    table_p = [570.375, 17.919, 3.375, 19.125, 3.375, 2.331, 3.375, 2.331, 3.375,
               19.125, 3.375, 17.919, 3.375]
    assert all(abs(a - b) < 5e-4 for a, b in zip(c, table_c))
    assert all(abs(a - b) < 5e-4 for a, b in zip(power, table_p))
    assert c[4] < c[3] and c[10] < c[9] and c[12] < c[11]           # contributions do not rise
    for g, share in [(1, 97.6), (2, 94.0), (4, 86.6)]:
        x = [((g * a) % 9) or 9 for a in r9]
        print("  class S=%4d  share %.1f%%" % (S(x), 100 * B(x) / S(x)))
        assert round(100 * B(x) / S(x), 1) == share

    print("== Section 5: maximality sweep and the spectral counterexamples, m <= 120 ==")
    wins, fails, B_fails, share_fails, not_unit = [], [], [], [], []
    for m in range(2, 121):
        cy = orbits(m)
        pr = orbit_of(cy, (0, 1))
        P = len(pr)
        if P % 2:
            continue
        rec = {}
        for c in cy:
            if len(c) == P:
                x = shown(c, m)
                rec[id(c)] = (S(x), B(x), c)
        smax = max(s for s, _, _ in rec.values())
        maxima = [v for v in rec.values() if v[0] == smax]
        if max(b for _, b, _ in maxima) < max(b for _, b, _ in rec.values()) - 1e-6:
            B_fails.append(m)
        s_pr, b_pr, _ = rec[id(pr)]
        if s_pr == smax:
            wins.append(m)
            continue
        fails.append(m)
        if max(b / s for s, b, _ in maxima) < max(b / s for s, b, _ in rec.values()) - 1e-9:
            share_fails.append(m)
        units = {(0, g) for g in range(1, m) if gcd(g, m) == 1}
        if not any(units & set(c) for _, _, c in maxima):
            not_unit.append(m)
        assert all(b > b_pr + 1e-6 for _, b, _ in maxima)           # strictly more B than (1,1)
        if m == 40:
            (s3, b3, _), = [v for v in maxima if (0, 3) in v[2]]
            s13, b13, _ = rec[id(orbit_of(cy, (0, 13)))]
            print("  m=40: (0,3) S=%d %.1f%% | (1,1) S=%d %.1f%% | (0,13) S=%d %.1f%%"
                  % (s3, 100 * b3 / s3, s_pr, 100 * b_pr / s_pr, s13, 100 * b13 / s13))
            assert (s3, s_pr, s13) == (58320, 51920, 51920)
            assert [round(100 * f, 1) for f in (b3 / s3, b_pr / s_pr, b13 / s13)] == [93.7, 88.8, 96.1]
        if m == 71:
            (s17, b17, _), = maxima
            assert (0, 17) in maxima[0][2]
            s29, b29, _ = rec[id(orbit_of(cy, (0, 29)))]
            print("  m=71: (0,17) S=%d B=%.1f | (0,29) S=%d B=%.1f" % (s17, b17, s29, b29))
            assert (s17, s29) == (242110, 240406)
            assert round(b17, 1) == 229974.2 and round(b29, 1) == 232216.4
        if m == 118:
            (s, b, _), = maxima
            print("  m=118: maximizer share %.2f%% vs (1,1) %.2f%%" % (100 * b / s, 100 * b_pr / s_pr))
            assert b / s < b_pr / s_pr - 1e-4
            assert (round(100 * b / s, 2), round(100 * b_pr / s_pr, 2)) == (93.94, 93.99)
    print("  even-Pisano moduli=%d  wins=%d  fails=%d" % (len(wins) + len(fails), len(wins), len(fails)))
    print("  failures:", fails)
    assert (len(wins), len(fails)) == (78, 40)
    assert fails == [15, 33, 40, 45, 49, 52, 53, 56, 59, 60, 64, 65, 67, 68, 71, 74, 80, 81,
                     83, 85, 86, 88, 90, 91, 92, 93, 97, 98, 103, 104, 105, 106, 107, 108,
                     109, 111, 115, 117, 118, 120]
    assert all(m in wins for m in range(3, 15)) and all(m in wins for m in (20, 25, 30, 35))
    decades = [sum(1 for m in fails if lo <= m < lo + 10) for lo in (80, 90, 100)]
    print("  failures per decade 80s/90s/100s:", decades)
    assert decades == [6, 6, 7]
    print("  maximizer misses the largest B at:", B_fails)
    assert B_fails == [70, 71, 75, 91, 102]
    print("  maximizer misses the largest share at %d of %d failures" % (len(share_fails), len(fails)))
    assert len(share_fails) == 24
    print("  failures with no unit-multiple maximizer (%d):" % len(not_unit), not_unit)
    assert len(not_unit) == 15

    print("== Section 6: S_max = m (P/2)^2 ==")
    def smax(m):
        cy = orbits(m)
        P = len(orbit_of(cy, (0, 1)))
        return P, max(S(shown(c, m)) for c in cy if len(c) == P)
    holds = []
    for m in range(2, 301):
        P, s = smax(m)
        if P % 2 == 0 and 4 * s == m * P * P:
            holds.append(m)
    print("  holds for even-period m <= 300 at:", holds)
    assert holds == [9, 27]
    for m, actual, predicted in [(3, 36, 48), (81, 857304, 944784), (243, 21677544, 25509168)]:
        P, s = smax(m)
        assert (s, m * (P // 2) ** 2) == (actual, predicted)
    P, s = smax(729)
    assert 4 * s != 729 * P * P
    print("  fails at 3, 81, 243 (values as stated) and at 729")

    print("== Section 6: the recurrence a_{n+1} = a_n + 2 a_{n-1}, odd m <= 101 ==")
    for m in range(3, 102, 2):
        cy = orbits(m, 1, 2)
        pr = orbit_of(cy, (0, 1))
        assert (1, 1) in pr
        same = [c for c in cy if len(c) == len(pr)]
        if m == 3:
            assert len(same) == 1
        else:
            assert S(shown(pr, m)) < max(S(shown(c, m)) for c in same)
    print("  (0,1) orbit is never the maximizer for 5 <= m <= 101; at m=3 it is the only orbit")

    print("\nAll Version 5 claims reproduced. OK.")


if __name__ == "__main__":
    main()
