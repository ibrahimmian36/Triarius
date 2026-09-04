#!/usr/bin/env python3
"""Session 102: the dense-misaligned observation as a theorem.

WHY THE EVENS DIE, worked out before coding. The evens admit no
coprime PAIR, so Theorem D cannot apply directly. But a class
(2m, a) with a EVEN kills only even x, and writing x = 2y its
condition becomes y = b (mod m): on the even cylinder the system
reduces to one whose moduli are ALL integers m, which does
contain a divergent coprime subfamily. Classes with a odd do the
same on the odd cylinder. Each cylinder's induced system is a
full system in disguise.

COROLLARY. If there is an M such that for every residue c mod M
the induced system at c contains a pairwise-coprime subfamily of
divergent drift, then dlog(A) exists and equals 0.
Proof: Theorem D gives that induced system density 0 AND
profinite measure 0, which is exactly the hypothesis of the
cylinder theorem; every cylinder contributes 0.

CHECKS.
 1. induced systems computed exactly (session-93 routine),
    elementwise;
 2. the misaligned dense families SATISFY the hypothesis - each
    cylinder's induced system has a coprime subfamily whose
    drift grows;
 3. THE ALIGNED families must FAIL it, since they survive. If
    they pass, the corollary is false. This check runs first.
"""
import math
import random
import sys
from math import gcd

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def induced(n, a, M, c):
    g = gcd(n, M)
    if (a - c) % g:
        return None
    nn = n // g
    if nn == 1:
        return (1, 0)
    inv = pow((M // g) % nn, -1, nn)
    return (nn, ((((a - c) // g) * inv) % nn))


def check_induced_exact(X=40000):
    ok, tested = True, 0
    # the guard fired on COUNT, not on a mismatch: only 6
    # non-empty pairs arose from five cases. Widen the sample.
    cases = [(n, a, M)
             for n in (6, 8, 10, 12, 15, 20, 21, 28)
             for a in (0, 1, 2, 5)
             for M in (2, 3, 4, 6) if a < n]
    for (n, a, M) in cases:
        for c in range(M):
            ind = induced(n, a, M, c)
            direct = {x for x in range(1, X)
                      if x % M == c % M and x % n == a % n}
            if ind is None:
                if direct: ok = False
                continue
            nn, f = ind
            # OFF-BY-ONE FIX: y must run to X//M + 1, else the
            # last few x below X are missed and every case
            # mismatches at the boundary (39996, 39997, ...).
            via = {c + M*y for y in range(0, X//M + 1)
                   if y % nn == f % nn and 0 < c + M*y < X}
            tested += 1
            if direct != via: ok = False
    ck("induced systems computed exactly", ok and tested > 8,
       "%d pairs" % tested)


def coprime_drift(mods, cap=4000):
    """greedy pairwise-coprime subfamily, drift"""
    chosen, d = [], 0.0
    for m in mods:
        if m < 2: continue
        if all(gcd(m, c) == 1 for c in chosen):
            chosen.append(m); d += 1.0/m
            if len(chosen) >= cap: break
    return d, len(chosen)


def hypothesis(cls, M, label):
    """does every cylinder's induced system have a coprime
    subfamily of growing drift?"""
    rows = []
    for c in range(M):
        ind = [induced(n, a, M, c) for n, a in cls]
        mods = [nn for x in ind if x is not None
                for nn in (x[0],) if nn > 1]
        d, k = coprime_drift(mods)
        rows.append((c, len(mods), k, d))
    print("   %-26s %s" % (label,
          "  ".join("c=%d: coprime drift %.3f (%d classes)"
                    % (c, d, k) for c, _, k, d in rows)))
    return all(d > 0.5 for _, _, _, d in rows)


def main():
    print("1. induced-system construction")
    check_induced_exact()
    print()
    rng = random.Random(7)
    LIM = 200000
    evens = [(n, rng.randrange(n)) for n in range(2, LIM, 2)]
    m3 = [(n, rng.randrange(n)) for n in range(3, LIM, 3)]
    ev_al = [(n, 0) for n in range(2, LIM, 2)]
    m3_al = [(n, 0) for n in range(3, LIM, 3)]

    print("2. ALIGNED families must FAIL the hypothesis "
          "(they survive) - checked first")
    a1 = hypothesis(ev_al, 2, "evens, residue 0")
    a2 = hypothesis(m3_al, 3, "multiples of 3, residue 0")
    ck("aligned families fail the hypothesis", not (a1 or a2),
       "evens pass=%s, mult3 pass=%s" % (a1, a2))
    print()
    print("3. MISALIGNED dense families should SATISFY it")
    b1 = hypothesis(evens, 2, "evens, random residues")
    b2 = hypothesis(m3, 3, "multiples of 3, random")
    ck("misaligned dense families satisfy the hypothesis",
       b1 and b2, "evens=%s, mult3=%s" % (b1, b2))
    print()
    if FAIL:
        print("FAILED: %s" % ", ".join(FAIL)); sys.exit(1)
    print("[verdict] the hypothesis separates exactly as the "
          "corollary requires: satisfied by the misaligned dense "
          "families that die, failed by the aligned ones that "
          "survive.")


if __name__ == "__main__":
    main()
