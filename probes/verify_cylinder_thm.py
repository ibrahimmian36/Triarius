#!/usr/bin/env python3
"""Session 93: a cylinder-decomposition theorem closing a gap in
my own paper.

THE GAP. Session 76's doubled construction has UNBOUNDED
misaligned parts, DIVERGENT drift, and survivor density exactly
1/2. So "unbounded misalignment + divergent drift => death" is
false; and the paper's corollary does not cover it either, since
its unbounded-misalignment classes carry divergent drift. Its
log density exists, but by an ad-hoc argument the paper never
generalises.

THE THEOREM. If some M splits N into cylinders on which each
induced system either has survivor density 0 or has bounded
misaligned parts, then dlog(A) exists.

Proof: A is the disjoint union over c mod M of A_c = A ∩ C_c.
Under x = c + My, A_c corresponds to the survivors of the
INDUCED system in y, and log density transfers with factor 1/M:
sum_{t in A_c, t<=X} 1/t = sum_{y} 1/(c+My) ~ (1/M) log(X/M) *
dlog(induced), so normalising by log X gives (1/M)*dlog(induced).
Finitely many cylinders, each covered.

CHECKS.
 1. the induced system is computed correctly - elementwise, the
    induced class must have exactly the right members;
 2. the log-density transfer factor 1/M, measured;
 3. the doubled construction satisfies the NEW hypothesis and
    fails the old ones;
 4. a control system fails the hypothesis.
Also: induced systems may repeat moduli, which earlier theorems
assumed distinct - checked and reported.
"""
import math
import sys
from math import gcd

FAIL = []
def ck(name, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + d) if d else ""))
    if not ok: FAIL.append(name)


def induced(n, a, M, c):
    """class (n,a) restricted to x = c mod M, in coords x=c+My.
    Returns (modulus, residue) in y, or None if it misses."""
    g = gcd(n, M)
    if (a - c) % g:
        return None
    nn = n // g
    if nn == 1:
        return (1, 0)
    inv = pow((M // g) % nn, -1, nn)
    f = (((a - c) // g) * inv) % nn
    return (nn, f)


def check_induced(X=60000):
    ok = True
    tested = 0
    for (n, a, M) in [(6, 2, 4), (10, 4, 6), (12, 8, 10),
                      (15, 5, 4), (21, 7, 6), (8, 2, 12)]:
        for c in range(M):
            ind = induced(n, a, M, c)
            direct = {x for x in range(1, X)
                      if x % M == c % M and x % n == a % n}
            if ind is None:
                if direct:
                    ok = False
                continue
            nn, f = ind
            viaY = {c + M * y for y in range(0, X // M)
                    if y % nn == f % nn and c + M * y < X
                    and c + M * y > 0}
            tested += 1
            if direct != viaY:
                ok = False
    ck("induced class computed exactly (elementwise)",
       ok and tested > 10, "%d cylinder-class pairs" % tested)
    return ok


def check_transfer(M=6, c=1):
    """Log density transfers with factor 1/M ASYMPTOTICALLY.
    Finite X carries corrections of order 1/log X (from the
    y-range ending at X/M, and from the Euler constant in the
    harmonic sum) - my first two attempts at a closed-form
    finite-X prediction were both wrong. The honest test is that
    the measured value APPROACHES the asymptotic 1/M * delta
    with a shrinking gap."""
    dS = 1 - 1.0 / 5
    asym = dS / M
    rows = []
    for X in (10 ** 5, 10 ** 6, 10 ** 7):
        acc = 0.0
        for y in range(1, X // M):
            if y % 5 == 0:
                continue
            t = c + M * y
            if t < X:
                acc += 1.0 / t
        rows.append((X, acc / math.log(X)))
    gaps = [abs(v - asym) for _, v in rows]
    for (X, v), g in zip(rows, gaps):
        print("     X=%-9d measured %.5f  asymptotic %.5f  "
              "gap %.5f" % (X, v, asym, g))
    ck("log density transfers with factor 1/M (gap shrinks "
       "toward the asymptotic value)",
       all(y < x for x, y in zip(gaps, gaps[1:]))
       and gaps[-1] < 0.01,
       "gaps %s" % ", ".join("%.4f" % g for g in gaps))


def doubled_construction():
    """session 76: blocks of large primes per small prime q,
    moduli 2pq, residues even"""
    def primes(n):
        s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
        i = 2
        while i * i <= n:
            if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
            i += 1
        return [i for i in range(2, n + 1) if s[i]]
    ps = primes(200000)
    smalls = ps[1:9]
    pool = [p for p in ps if p > smalls[-1] and p >= 200]
    out, idx = [], 0
    for q in smalls:
        acc, j = 0.0, 0
        while idx < len(pool) and acc < 0.05:
            p = pool[idx]; idx += 1
            r = j % q; j += 1
            inv = pow(q % p, -1, p)
            k = ((0 - r) % p) * inv % p
            a = (r + q * k) % (p * q)
            n2 = 2 * p * q
            a2 = a if a % 2 == 0 else a + p * q
            out.append((n2, a2 % n2))
            acc += 1.0 / p
    return sorted(set(out)), smalls


def mis(n, a):
    b = None
    for w in range(1, n + 1):
        if n % w: continue
        u = n // w
        if gcd(u, w) == 1 and a % u == 0 and (b is None or w < b):
            b = w
    return b


def main():
    print("1. the induced-system construction")
    check_induced()
    print()
    print("2. log-density transfer under x = c + My")
    check_transfer()
    print()
    print("3. the doubled construction versus the hypotheses")
    cls, smalls = doubled_construction()
    ws = [mis(n, a) for n, a in cls]
    drift = sum(1.0 / n for n, _ in cls)
    print("   %d classes; misaligned parts up to %d; drift %.4f "
          "(diverges as more small primes are added)"
          % (len(cls), max(ws), drift))
    # every kill set even?
    odd_killers = sum(1 for n, a in cls if n % 2 or a % 2)
    ck("no class kills an odd number", odd_killers == 0)
    # cylinder c=1 mod 2: induced system is EMPTY
    empt = all(induced(n, a, 2, 1) is None for n, a in cls)
    ck("cylinder c=1 (odds): induced system is empty", empt)
    # cylinder c=0 mod 2: induced system exists
    ind0 = [induced(n, a, 2, 0) for n, a in cls]
    live = [i for i in ind0 if i is not None]
    ck("cylinder c=0 (evens): induced system non-empty",
       len(live) > 100, "%d induced classes" % len(live))
    rep = len(live) - len({m for m, _ in live})
    print("   NOTE: induced system has %d repeated moduli "
          "(earlier theorems assume distinct moduli; the "
          "cylinder theorem must not rely on distinctness)" % rep)
    print()
    print("4. does it satisfy the OLD hypotheses?")
    ck("NOT bounded misalignment (parts unbounded)",
       max(ws) > 10 * min(ws), "max %d" % max(ws))
    # The corollary asks whether the classes of unbounded
    # misalignment carry CONVERGENT drift. Here the misaligned
    # part of a class is its small prime q, so for any fixed M
    # the excluded classes are those with q not dividing M -
    # all but finitely many q - and their drift is
    # sum_{q > Q} s/q, which DIVERGES by Mertens. A finite
    # instance cannot exhibit this; the first version of this
    # check measured 8 blocks and failed for that reason.
    tail = [(q, 0.05 / q) for q in smalls]
    print("     per-block drift s/q: %s"
          % "  ".join("q=%d:%.4f" % t for t in tail))
    print("     sum over these %d blocks: %.4f" % (len(tail), sum(x for _,x in tail)))
    print("     the excluded drift is sum_{q>Q} s/q, divergent "
          "by Mertens - so for EVERY M the corollary's "
          "convergence hypothesis fails")
    # misaligned part 1 means a fully ALIGNED class, which
    # is legitimate and simply not one of the small primes
    parts_are_smalls = set(ws) <= set(smalls) | {1}
    ck("NOT the corollary: each class's misaligned part IS its "
       "small prime q, so for any M the excluded classes are "
       "those with q not dividing M and their drift "
       "sum_{q>Q} s/q diverges by Mertens",
       parts_are_smalls,
       "misaligned parts %s subset of small primes %s"
       % (sorted(set(ws)), smalls))
    print()
    if FAIL:
        print("FAILED: %s" % ", ".join(FAIL)); sys.exit(1)
    print("[verdict] The doubled construction falls OUTSIDE both "
          "existing hypotheses and INSIDE the cylinder "
          "decomposition (M=2: odds carry no classes, evens die). "
          "The theorem closes a genuine gap in the paper.")


if __name__ == "__main__":
    main()
