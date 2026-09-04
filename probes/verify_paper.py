#!/usr/bin/env python3
"""Session 82: independent verification of every checkable claim
in papers/erdos25/main.tex.

Written FRESH rather than by importing the probes that produced
the claims: re-running the same code proves only determinism,
whereas independent reimplementation can actually catch a wrong
claim. Where a quantity has both a closed form and a direct
computation, both are computed and compared.

Registered: every claim either reproduces, or the PAPER is
edited to match reality - never the reverse. A claim that cannot
be reproduced at all is cut from the paper, not softened.
"""

import sys
from fractions import Fraction
from itertools import combinations, product
from math import gcd, log

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  " + detail) if detail else ""))
    if not ok:
        FAIL.append(name)
    return ok


# ---------- Lemma A: incompatible drift <= 1, sharp ----------
def lemma_A():
    U = [(m, b) for m in range(2, 15) for b in range(m)]
    n = len(U)
    inc = [[True] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            (m1, b1), (m2, b2) = U[i], U[j]
            inc[i][j] = ((b1 - b2) % gcd(m1, m2) != 0)
    order = sorted(range(n), key=lambda i: U[i][0])
    best = [Fraction(0)]
    w = [Fraction(1, U[i][0]) for i in order]
    suf = [Fraction(0)] * (n + 1)
    for i in range(n - 1, -1, -1):
        suf[i] = suf[i + 1] + w[i]

    def rec(i, chosen, acc):
        if acc > best[0]:
            best[0] = acc
        if i >= n or acc + suf[i] <= best[0]:
            return
        v = order[i]
        if all(inc[v][u] for u in chosen):
            chosen.append(v)
            rec(i + 1, chosen, acc + w[i])
            chosen.pop()
        rec(i + 1, chosen, acc)

    rec(0, [], Fraction(0))
    ok1 = check("Lemma A: exact max over 104 classes (mod<=14) "
                "equals 1", best[0] == 1 and n == 104,
                "max=%s, classes=%d" % (best[0], n))
    # the extremal family partitions N
    fam = [(2, 0)] + [(12, r) for r in (1, 3, 5, 7, 9, 11)]
    drift = sum(Fraction(1, m) for m, _ in fam)
    cover = all(any(x % m == b % m for m, b in fam)
                for x in range(1, 3000))
    pairwise = all((b1 - b2) % gcd(m1, m2) != 0
                   for (m1, b1), (m2, b2)
                   in combinations(fam, 2))
    ok2 = check("Lemma A: extremal family is incompatible, has "
                "drift 1, and partitions N",
                drift == 1 and cover and pairwise,
                "drift=%s covers=%s incompatible=%s"
                % (drift, cover, pairwise))
    return ok1 and ok2


# ---------- Theorem H and Theorem D ----------
def density(classes, X):
    k = bytearray(X)
    for n, a in classes:
        s = a % n
        if s == 0:
            s = n
        while s < n:
            s += n
        for i in range(s, X, n):
            k[i] = 1
    lo = X // 2
    return 1.0 - sum(k[lo:X]) / float(X - lo)


def thm_H_D(X=600000):
    fam = [(2 ** j, 2 ** (j - 1)) for j in range(1, 15)]
    inc = all((b1 - b2) % gcd(m1, m2) != 0
              for (m1, b1), (m2, b2) in combinations(fam, 2))
    pred = 1.0 - sum(1.0 / n for n, _ in fam)
    got = density(fam, X)
    okH = check("Theorem H: incompatible => density = 1 - sum 1/n",
                inc and abs(pred - got) < 0.002,
                "incompatible=%s predicted %.5f measured %.5f"
                % (inc, pred, got))
    ps = [p for p in range(2, 200)
          if all(p % q for q in range(2, int(p ** .5) + 1))]
    cop = [(p, 1 % p) for p in ps]
    predD = 1.0
    for n, _ in cop:
        predD *= (1 - 1.0 / n)
    gotD = density(cop, X)
    okD = check("Theorem D: coprime => density = prod(1-1/n)",
                abs(predD - gotD) < 0.005,
                "predicted %.5f measured %.5f" % (predD, gotD))
    mut = [(2 * p, 2 % (2 * p)) for p in ps[1:14]]
    pm = 1.0
    for n, _ in mut:
        pm *= (1 - 1.0 / n)
    gm = density(mut, X)
    okM = check("  mutation: non-coprime family must DEVIATE",
                abs(pm - gm) > 0.02,
                "deviation %.4f" % abs(pm - gm))
    return okH and okD and okM


# ---------- Proposition: the D-E break ----------
def de_break(X=100):
    # PAPER CORRECTED: the counterexample requires the FINITE
    # head. In the full system (n,1) every x >= 3 is killed by
    # the class n = x-1, so Lemma 1 holds with equality and
    # there is no violation. Session 79 used moduli 2..39.
    cls = [(n, 1 % n) for n in range(2, 40)]
    killed = set()
    for n, a in cls:
        s = a % n
        if s == 0:
            s = n
        while s < n:
            s += n
        for x in range(s, 3000, n):
            killed.add(x)
    ok1 = check("Prop: 31 is killed and 62 is not (closure fails)",
                31 in killed and 62 not in killed,
                "31 killed=%s, 62 killed=%s"
                % (31 in killed, 62 in killed))
    # Lemma 1 at n = 42
    def Lam(k):
        for p in range(2, k + 1):
            if k % p == 0:
                q, e = p, 0
                m = k
                while m % p == 0:
                    m //= p
                    e += 1
                return log(p) if m == 1 else 0.0
        return 0.0
    divs = [d for d in range(1, 43) if 42 % d == 0]
    rhs = sum(Lam(42 // d) for d in divs if d in killed)
    lhs = (1.0 if 42 in killed else 0.0) * log(42)
    expect = log(2) + log(3) + log(7)
    ok2 = check("Prop: Lemma 1 fails at n=42 with lhs 0 and rhs "
                "log2+log3+log7", lhs == 0.0
                and abs(rhs - expect) < 1e-9,
                "lhs=%.3f rhs=%.3f expected %.3f"
                % (lhs, rhs, expect))
    return ok1 and ok2


# ---------- Theorem: cylinder reduction ----------
def cylinder(X=40000):
    def mis(n, a):
        best = None
        for w in range(1, n + 1):
            if n % w:
                continue
            u = n // w
            if gcd(u, w) == 1 and a % u == 0:
                if best is None or w < best:
                    best = w
        return best
    cls = [(6, 2), (10, 4), (14, 6), (15, 5), (21, 7)]
    ws = [mis(n, a) for n, a in cls]
    M = 1
    for w in ws:
        M = M * w // gcd(M, w)
    bad = tested = 0
    for c in range(M):
        for (n, a), w in zip(cls, ws):
            u = n // w
            S = set(x for x in range(1, X)
                    if x % M == c and x % n == a % n)
            T = set(x for x in range(1, X)
                    if x % M == c and x % u == 0)
            if not S:
                continue
            tested += 1
            if S != T:
                bad += 1
    return check("Theorem (bounded misalignment): cylinder "
                 "reduction exact, misaligned parts %s, M=%d"
                 % (ws, M), bad == 0 and tested > 20,
                 "%d/%d pairs exact" % (tested - bad, tested))


# ---------- Heilbronn-Rohrbach failure: RE-DERIVE ----------
def hr(X=60000):
    """the paper cites '42.4% of 55,770 systems'. The sweep's
    parameter space was never recorded, so re-derive it here
    over a STATED space and correct the paper to match."""
    tri = list(combinations(range(2, 13), 3))
    total = viol = 0
    worst = None
    for (a, b, c) in tri:
        for ra in range(a):
            for rb in range(b):
                for rc in range(c):
                    cls = [(a, ra), (b, rb), (c, rc)]
                    total += 1
                    bound = 1.0
                    for n, _ in cls:
                        bound *= (1 - 1.0 / n)
                    d = density(cls, X)
                    if d < bound - 1e-9:
                        viol += 1
                        if worst is None or d - bound < worst[0]:
                            worst = (d - bound, cls, d, bound)
    pct = 100.0 * viol / total
    print("  [INFO] HR sweep over all 3-class systems with "
          "moduli from {2..12}: %d systems, %d violate (%.1f%%)"
          % (total, viol, pct))
    print("         worst: %s density %.4f vs bound %.6f"
          % (worst[1], worst[2], worst[3]))
    named = [(2, 0), (4, 1), (8, 3)]
    dn = density(named, X)
    bn = 1.0
    for n, _ in named:
        bn *= (1 - 1.0 / n)
    ok = check("HR: named counterexample (2,0),(4,1),(8,3)",
               abs(dn - 0.125) < 0.002 and abs(bn - 0.328125) < 1e-9,
               "density %.4f vs bound %.6f" % (dn, bn))
    return ok, total, pct


def main():
    print("Independent verification of papers/erdos25/main.tex")
    print()
    lemma_A()
    thm_H_D()
    de_break()
    cylinder()
    _, hr_total, hr_pct = hr()
    print()
    if FAIL:
        print("FAILED CLAIMS: %s" % ", ".join(FAIL))
        print("-> the PAPER must be corrected to match reality.")
        sys.exit(1)
    print("All reproduced claims PASS.")
    print("HR figures now in the paper: %d systems, %.1f%%"
          % (hr_total, hr_pct))
    # confirm the full system really does kill everything,
    # which is WHY the finite head is required
    full_killed = all(any(x % n == 1 % n and x >= n
                          for n in range(2, x + 1))
                      for x in range(3, 400))
    print("cross-check: in the FULL system (n,1) every x>=3 is "
          "killed -> %s (this is why the paper now states the "
          "finite head)" % full_killed)


if __name__ == "__main__":
    main()
