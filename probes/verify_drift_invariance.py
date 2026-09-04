#!/usr/bin/env python3
"""Session 95: does cylinder refinement starve cylinders of
drift?  A closed-form identity says no.

DERIVATION. Fix M. A class (n,a) meets the cylinder c mod M iff
g | a-c where g = gcd(n,M), so it meets exactly M/g of the M
cylinders. On such a cylinder it becomes, in the coordinate
x = c+My, a class of modulus n/g, contributing induced drift
g/n there. Its total contribution over all cylinders is
therefore (M/g)*(g/n) = M/n, independent of M and of g. Summing
over classes, the total induced drift over all cylinders is
M * sum_i 1/n_i, so

    average induced drift per cylinder = sum_i 1/n_i,

exactly the ORIGINAL drift, for every M.

CONSEQUENCE. Refinement preserves drift on average; it cannot
resolve cylinders by starving them, and a divergent drift stays
divergent on the typical cylinder. This closes off the natural
hope behind iterating the cylinder decomposition.

CHECK: compute the induced drift on every cylinder directly and
compare the average with the original drift, exactly, for a
range of M and several systems. CONTROL: a system of known
drift. The identity must hold for EVERY M tested, not one.
"""
import sys
from fractions import Fraction
from math import gcd

FAIL = []
def ck(name, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + d) if d else ""))
    if not ok: FAIL.append(name)


def induced_drift_sum(cls, M):
    """sum over all cylinders of the induced drift, exactly"""
    tot = Fraction(0)
    for n, a in cls:
        g = gcd(n, M)
        for c in range(M):
            if (a - c) % g == 0:
                tot += Fraction(g, n)
    return tot


def systems():
    S = {}
    S["aligned 2..60"] = [(n, 0) for n in range(2, 60)]
    S["shifted 2..60"] = [(n, 1 % n) for n in range(2, 60)]
    import random
    rng = random.Random(5)
    S["random 2..60"] = [(n, rng.randrange(n))
                         for n in range(2, 60)]
    S["even moduli"] = [(2 * k, k % (2 * k))
                        for k in range(1, 40)]
    S["prime moduli"] = [(p, 3 % p) for p in
                         (2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                          31, 37, 41, 43, 47)]
    return S


def main():
    print("average induced drift per cylinder vs original drift")
    print()
    Ms = [1, 2, 3, 4, 6, 8, 12, 30, 60, 210]
    allok = True
    for name, cls in sorted(systems().items()):
        orig = sum(Fraction(1, n) for n, _ in cls)
        bad = []
        for M in Ms:
            avg = induced_drift_sum(cls, M) / M
            if avg != orig:
                bad.append((M, avg))
        print("  %-16s original drift %.6f   %s"
              % (name, float(orig),
                 "identity holds for all %d moduli" % len(Ms)
                 if not bad else "FAILS at %s" % bad[:2]))
        allok &= not bad
    print()
    ck("average induced drift equals the original drift, for "
       "every tested M and every system", allok,
       "%d systems x %d moduli" % (len(systems()), len(Ms)))

    # control: a system whose drift is known by hand
    ctrl = [(2, 0), (3, 1), (5, 2)]
    known = Fraction(1, 2) + Fraction(1, 3) + Fraction(1, 5)
    ok = all(induced_drift_sum(ctrl, M) / M == known
             for M in Ms)
    ck("control with hand-computed drift 1/2+1/3+1/5", ok,
       "= %s" % known)

    print()
    print("  Consequence: refinement does NOT dilute drift.  A")
    print("  cylinder decomposition cannot resolve cylinders by")
    print("  starving them, and divergent drift remains")
    print("  divergent on the typical cylinder.")
    print()
    if FAIL:
        print("FAILED: %s" % ", ".join(FAIL)); sys.exit(1)
    print("[verdict] identity confirmed exactly; the starvation "
          "mechanism behind iterated cylinder decomposition is "
          "closed off.")


if __name__ == "__main__":
    main()
