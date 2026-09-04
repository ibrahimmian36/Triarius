#!/usr/bin/env python3
"""Session 73: PROVE Conjecture C on structured subfamilies.

RESEARCH FINDING THAT SHAPES THIS PROBE. C is FALSE as an
abstract weighted-graph statement: take m disjoint cliques, each
k vertices of weight B/k. Every clique has weight B, every
independent set (one vertex per clique) has weight m*B/k <= 1,
yet the total is m*B ~ k, unbounded. So no proof of C can come
from Lemma A alone - it MUST use arithmetic. The arithmetic that
kills that construction is distinctness of moduli: the block
[M,2M) holds at most M distinct moduli and so carries at most
log 2 of drift, while the construction needs k^2/B distinct
moduli all of size ~k/B. Hence: prove C on arithmetically
structured subfamilies.

LEMMA C1 (coprime). Pairwise-coprime moduli are pairwise
COMPATIBLE: gcd = 1 makes the congruence condition vacuous. So a
pairwise-coprime subfamily with divergent drift IS a divergent
channel, and C holds for it.

LEMMA C3 (constant-gcd pigeonhole). Let F have gcd(n,n') = d for
every distinct pair. Partition F by residue mod d. Within one
part, any two classes agree mod d = gcd, so they are compatible:
each part is a CHANNEL, and there are exactly d parts. If F has
divergent drift then some part has divergent drift, since a
finite sum of convergent series converges. So C holds for F.
C1 is C3 at d = 1.

This proves C whenever divergence is carried by a common-cofactor
family - which is exactly what the session-72 Q-families were,
and explains the (1/Q) log log K they produced.

CHECKS. (1) C1 exactly. (2) C3's parts are genuinely channels,
verified pairwise, with the MUTATION that classes differing mod d
must be INCOMPATIBLE - else the partition is doing no work.
(3) the pigeonhole is quantitative: the best part must carry at
least 1/d of the family's drift. All exact integer arithmetic.
Pod-only, stdlib only.
"""

import sys
import time
from math import gcd


def compat(c1, c2):
    (m1, b1), (m2, b2) = c1, c2
    return (b1 - b2) % gcd(m1, m2) == 0


def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    i = 2
    while i * i <= n:
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
        i += 1
    return [i for i in range(2, n + 1) if s[i]]


def check_C1():
    """pairwise coprime => pairwise compatible, for ALL residues"""
    ps = primes_upto(120)
    bad = 0
    tested = 0
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            for b1 in range(min(ps[i], 7)):
                for b2 in range(min(ps[j], 7)):
                    tested += 1
                    if not compat((ps[i], b1), (ps[j], b2)):
                        bad += 1
    print("[1] LEMMA C1: coprime moduli are compatible for EVERY "
          "residue pair -- %d pairs tested, %d failures (must be "
          "0)" % (tested, bad))
    return bad == 0 and tested > 1000


def check_C3(ds=(2, 4, 6, 12, 30), nprimes=60):
    """constant-gcd families: residue-mod-d parts are channels,
    and classes in DIFFERENT parts must be incompatible"""
    ps = primes_upto(400)
    ok_parts = True
    mut_fired = 0
    quant_ok = True
    for d in ds:
        # moduli d*p for primes p not dividing d -> pairwise
        # gcd is exactly d
        fam = [d * p for p in ps if d % p and p > d][:nprimes]
        gset = set()
        for i in range(len(fam)):
            for j in range(i + 1, len(fam)):
                gset.add(gcd(fam[i], fam[j]))
        if gset != {d}:
            print("    [setup] d=%d: pairwise gcds are %s, not "
                  "{%d} -- family invalid" % (d, sorted(gset), d))
            return False
        # assign residues: part index r means residue = r mod d
        parts = {}
        for idx, n in enumerate(fam):
            r = idx % d
            b = r % d
            # lift r to a residue mod n congruent to r mod d
            while b % d != r:
                b += 1
            parts.setdefault(r, []).append((n, b))
        # every part is a channel
        for r, cls in parts.items():
            for i in range(len(cls)):
                for j in range(i + 1, len(cls)):
                    if not compat(cls[i], cls[j]):
                        ok_parts = False
        # MUTATION: different parts must be INCOMPATIBLE
        keys = sorted(parts)
        for x in range(len(keys)):
            for y in range(x + 1, len(keys)):
                a, b = parts[keys[x]][0], parts[keys[y]][0]
                if not compat(a, b):
                    mut_fired += 1
        # quantitative pigeonhole: best part >= (1/d) of total
        tot = sum(1.0 / n for n in fam)
        best = max(sum(1.0 / n for n, _ in cls)
                   for cls in parts.values())
        if best < tot / d - 1e-12:
            quant_ok = False
    print("[2] LEMMA C3: every residue-mod-d part is a channel "
          "(pairwise compatible): %s" % ("PASS" if ok_parts
                                         else "FAIL"))
    print("    MUTATION: classes in DIFFERENT parts are "
          "incompatible in %d checked cross-pairs (need >= 1, "
          "else the partition does nothing)" % mut_fired)
    print("[3] quantitative pigeonhole, best part >= (1/d) of "
          "the family drift: %s" % ("PASS" if quant_ok
                                    else "FAIL"))
    return ok_parts and mut_fired >= 1 and quant_ok


def check_finite_support():
    """divergent drift forces INFINITE prime support: over a
    finite prime set S the drift is bounded by the Euler product
    prod_{p in S} (1 - 1/p)^-1 < infinity"""
    S = [2, 3, 5, 7, 11]
    lim = 10 ** 7
    tot, stack = 0.0, [1]
    seen = set()
    while stack:
        x = stack.pop()
        if x in seen or x > lim:
            continue
        seen.add(x)
        if x >= 2:
            tot += 1.0 / x
        for p in S:
            if x * p <= lim:
                stack.append(x * p)
    euler = 1.0
    for p in S:
        euler *= 1.0 / (1.0 - 1.0 / p)
    print("[4] finite prime support {2,3,5,7,11}: drift over all "
          "such moduli <= %.4f, Euler bound %.4f -- CONVERGENT, "
          "so divergent drift forces infinite prime support"
          % (tot, euler - 1.0))
    return tot <= euler


def main():
    t0 = time.time()
    r = [check_C1(), check_C3(), check_finite_support()]
    print()
    if all(r):
        print("[C1 and C3 VERIFIED] Conjecture C is PROVED for "
              "every system whose divergence is carried by a "
              "constant-gcd (common-cofactor) family, coprime "
              "families included as d = 1. The residual case is "
              "divergence spread over infinitely many cofactors, "
              "each individually convergent.")
    else:
        print("[FAILED] %s" % r)
        sys.exit(1)
    print("[done] %ds" % (time.time() - t0))


if __name__ == "__main__":
    main()
