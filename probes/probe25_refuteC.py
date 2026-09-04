#!/usr/bin/env python3
"""Session 74: an attempted REFUTATION of Conjecture C.

GEOMETRIC REFRAMING (session 74 research). A class (n,a) is the
ball a + n*Zhat in the profinite integers, of Haar measure 1/n.
Two classes are COMPATIBLE exactly when their balls intersect,
and since p-adic balls are nested-or-disjoint, a pairwise
compatible family has a COMMON POINT. So:

    a channel = a point c in Zhat lying in all its balls.

C says: sum 1/n_i = infinity implies some point lies in a
divergent-measure subfamily.

TWO IMMEDIATE CONSEQUENCES.
 * For PRIME moduli, CRT plus compactness gives a c with
   c = a_p mod p for EVERY p at once, so the channel holds all
   prime classes and has drift sum 1/p = infinity. This explains
   session 71's measured 0.998 log log K exactly: the greedy
   channel IS the prime sum.
 * integral of S(c) = sum 1/n_i^2, which CONVERGES. So almost
   every c lies in only finitely much mass, and C needs a
   measure-ZERO point. No averaging or probabilistic argument
   can prove C.

THE CONSTRUCTION. Since the standard family satisfies C only
through its primes, remove them. Partition the large primes into
consecutive blocks P_q of harmonic weight ~s, one block per
small prime q, and take moduli n = p*q for p in P_q. Then:
  - different blocks: gcd(pq, p'q') = 1, so ALWAYS compatible;
  - same block: gcd(pq, p'q) = q exactly, so compatible iff the
    residues agree mod q.
Choose residues spread uniformly mod q. The maximum channel then
DECOMPOSES exactly: pick, from each block, the heaviest residue
group, and take the union (legal, since blocks are pairwise
coprime). So

  total drift   = sum_q s/q          ~ s * log log Q -> INFINITY
  max channel   = sum_q (best group) ~ sum_q s/q^2   -> FINITE

If that holds, C is FALSE. The maximum is computed EXACTLY from
the decomposition, not by greedy; greedy is run only as an
independent lower bound and must never exceed the exact value.

MUTATION (must fire): aligning the residues within each block
collapses the block into ONE channel, making max channel = total
drift. If that does not happen, the adversarial residue choice
is not load-bearing and the construction proves nothing.
Pod-only, stdlib only, exact integer arithmetic.
"""

import argparse
import math
import random
import sys
import time
from math import gcd


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= n:
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
        i += 1
    return [i for i in range(2, n + 1) if s[i]]


def crt(r1, m1, r2, m2):
    """solve x = r1 mod m1, x = r2 mod m2 with gcd(m1,m2)=1"""
    inv = pow(m1 % m2, -1, m2)
    k = ((r2 - r1) % m2) * inv % m2
    return r1 + m1 * k


def build(nblocks, s_target, plimit, mode, pstart=200):
    """blocks of large primes, weight ~s_target each, one per
    small prime q; moduli p*q; residues spread mod q (mode
    'spread') or all equal (mode 'aligned', the mutation)"""
    ps = primes_upto(plimit)
    smalls = ps[:nblocks]
    big_start = 0
    while big_start < len(ps) and ps[big_start] <= smalls[-1]:
        big_start += 1
    pool = [p for p in ps[big_start:] if p >= pstart]
    blocks, idx = [], 0
    for q in smalls:
        acc, blk = 0.0, []
        while idx < len(pool) and acc < s_target:
            p = pool[idx]
            idx += 1
            if p == q:
                continue
            blk.append(p)
            acc += 1.0 / p
        if acc < s_target * 0.98:
            return None, None, None   # ran out of primes
        blocks.append((q, blk, acc))
    classes, per_block = [], []
    for bi, (q, blk, acc) in enumerate(blocks):
        groups = {}
        for j, p in enumerate(blk):
            r = 0 if mode == "aligned" else j % q
            a = crt(r % q, q, 0, p)          # a = r mod q, 0 mod p
            n = p * q
            assert a % q == r % q and a % p == 0
            classes.append((n, a % n))
            groups.setdefault(r % q, 0.0)
            groups[r % q] += 1.0 / n
        per_block.append((q, acc, groups))
    return classes, per_block, blocks


def verify_structure(classes, blocks):
    """exact checks the whole argument rests on"""
    ns = [n for n, _ in classes]
    if len(set(ns)) != len(ns):
        return "moduli not distinct"
    qof = {}
    for q, blk, _ in blocks:
        for p in blk:
            qof[p * q] = q
    cross_checked = same_checked = 0
    rng = random.Random(5)
    sample = classes if len(classes) <= 400 else \
        rng.sample(classes, 400)
    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            (n1, a1), (n2, a2) = sample[i], sample[j]
            q1, q2 = qof[n1], qof[n2]
            g = gcd(n1, n2)
            if q1 != q2:
                if g != 1:
                    return "cross-block gcd %d, expected 1" % g
                cross_checked += 1
            else:
                if g != q1:
                    return "same-block gcd %d, expected %d" % (g, q1)
                ok = (a1 - a2) % q1 == 0
                if ok != ((a1 - a2) % g == 0):
                    return "compatibility mismatch"
                same_checked += 1
    if cross_checked < 10 or same_checked < 10:
        return "checks vacuous (cross %d, same %d)" % (
            cross_checked, same_checked)
    return None


def exact_max_channel(per_block):
    """the decomposition: union over blocks of the heaviest
    residue group. Legal because distinct blocks are coprime,
    hence unconditionally compatible."""
    return sum(max(g.values()) for _, _, g in per_block)


def greedy_lower_bound(classes, restarts, seed):
    rng = random.Random(seed)
    best = 0.0
    cl = list(classes)
    for t in range(restarts):
        order = cl[:] if t == 0 else cl[:]
        if t:
            rng.shuffle(order)
        chosen, drift = [], 0.0
        for n, a in order:
            if all((a - b) % gcd(n, m) == 0 for m, b in chosen):
                chosen.append((n, a))
                drift += 1.0 / n
        best = max(best, drift)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plimit", type=int, default=4000000)
    ap.add_argument("--s", type=float, default=0.05)
    ap.add_argument("--pstart", type=int, default=200)
    ap.add_argument("--restarts", type=int, default=8)
    args = ap.parse_args()
    t0 = time.time()
    print("[construction] blocks of large primes of harmonic "
          "weight ~%.2f, one per small prime; moduli p*q"
          % args.s)
    print()
    rows = []
    for nb in (3, 5, 8, 12, 16, 20, 25):
        classes, per_block, blocks = build(nb, args.s,
                                           args.plimit, "spread",
                                           args.pstart)
        if classes is None:
            print("   %d blocks: not enough primes below %d "
                  "-- stopping" % (nb, args.plimit))
            break
        err = verify_structure(classes, blocks)
        if err:
            print("[STRUCTURE FAIL] %s" % err)
            sys.exit(1)
        total = sum(1.0 / n for n, _ in classes)
        exact = exact_max_channel(per_block)
        rows.append((nb, total, exact, len(classes)))
        print("   %d blocks | %6d classes | total drift %.4f | "
              "EXACT max channel %.4f | ratio %.3f"
              % (nb, len(classes), total, exact, exact / total))
    print()
    # independent cross-check on the largest feasible instance
    nb = rows[-1][0]
    classes, per_block, blocks = build(nb, args.s, args.plimit,
                                       "spread", args.pstart)
    sub = classes if len(classes) <= 1200 else \
        random.Random(1).sample(classes, 1200)
    sub_pb = None
    g = greedy_lower_bound(sub, args.restarts, 0)
    exact_sub = exact_max_channel(
        [(q, acc, {r: v for r, v in gr.items()})
         for q, acc, gr in per_block])
    print("[cross-check] restart-greedy on a %d-class subsample "
          "gives %.4f; the exact maximum for the full system is "
          "%.4f. Greedy must NOT exceed it: %s"
          % (len(sub), g, exact_sub,
             "OK" if g <= exact_sub + 1e-9 else "VIOLATED"))
    if g > exact_sub + 1e-9:
        print("[FAIL] greedy beat the claimed exact maximum, so "
              "the decomposition is wrong")
        sys.exit(1)
    print()
    # MUTATION
    cls_a, pb_a, blk_a = build(nb, args.s, args.plimit,
                               "aligned", args.pstart)
    tot_a = sum(1.0 / n for n, _ in cls_a)
    ex_a = exact_max_channel(pb_a)
    print("[mutation] with residues ALIGNED inside each block, "
          "max channel %.4f vs total %.4f (ratio %.3f) -- must be "
          "~1, else the spread residues are not load-bearing"
          % (ex_a, tot_a, ex_a / tot_a))
    if ex_a / tot_a < 0.99:
        print("[FAIL] mutation did not fire")
        sys.exit(1)
    print()
    ratios = [r[2] / r[1] for r in rows]
    falling = all(y < x for x, y in zip(ratios, ratios[1:]))
    print("   channel/total ratio by block count: %s"
          % "  ".join("%.3f" % x for x in ratios))
    print("   monotonically falling: %s" % falling)
    print()
    print("ANALYTIC EXTRAPOLATION (the part numerics cannot "
          "reach). Total drift = sum_q s/q over small primes q, "
          "which by Mertens DIVERGES like s*log log Q. Max "
          "channel = sum_q (best group) <= sum_q s/q^2, which "
          "CONVERGES (bounded by s*0.4523, the prime zeta "
          "P(2)). The two sums are computed over the same q, so "
          "the gap is unbounded.")
    print()
    if falling:
        print("[verdict] CONJECTURE C: REFUTATION CANDIDATE. The "
              "construction has divergent total drift with a "
              "maximum channel bounded by s*P(2), computed "
              "EXACTLY from the block decomposition, cross-"
              "checked against restart-greedy, with the residue "
              "mutation firing. CANDIDATE ONLY: this needs exact "
              "re-derivation on paper and the hostile gauntlet "
              "before any claim.")
    else:
        print("[verdict] construction does NOT separate: the "
              "ratio is not falling, so C survives this attack.")
    print("[done] %ds" % (time.time() - t0))


if __name__ == "__main__":
    main()
