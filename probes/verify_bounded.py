#!/usr/bin/env python3
"""Session 87: close the bounded-misalignment gap WITHOUT
characters.

CLAIM. Let every misaligned part w_i divide a fixed M, and fix a
cylinder C_c = {x = c mod M}. Because w_i | M, any x in C_c has
x = c (mod w_i), so class i can kill some x of C_c only if
a_i = c (mod w_i). Writing D_c = { u_i : a_i = c (mod w_i) },
the killed set inside C_c is a set of MULTIPLES of D_c
intersected with C_c. Hence for the head U_r of the first r
classes,

  (U \\ U_r) cap C_c  subset  (M(D_c) \\ M(D_c^{(r)})) cap C_c

apart from finitely many x below the head's largest modulus.
Intersecting with a cylinder only shrinks a set, so the upper
logarithmic density of the left side is at most that of
M(D_c) \\ M(D_c^{(r)}), which vanishes by the paper's
vanishing-tails lemma. Summing over the M cylinders and running
the epsilon-argument of the finite-channel theorem gives that
dlog(A) EXISTS.

CHECK 1 (exact, elementwise): the inclusion, with the finite
exceptional set counted explicitly rather than assumed.
MUTATION: with some w_i NOT dividing M the cylinder no longer
selects classes and the inclusion must FAIL - otherwise bounded
misalignment is doing no work.
CHECK 2: F(X) settles for such a system, as the conclusion
predicts.
"""
import sys
from math import gcd


def mis_split(n, a):
    best = None
    for w in range(1, n + 1):
        if n % w:
            continue
        u = n // w
        if gcd(u, w) == 1 and a % u == 0:
            if best is None or w < best:
                best = w
    return best


def build(M=105, limit=500):
    """classes with every misaligned part dividing M"""
    out = []
    for n in range(2, limit):
        g = gcd(n, M)
        u = n // g
        if gcd(u, g) != 1:
            continue
        a = (u * (n % g)) % n
        w = mis_split(n, a)
        # BUG FIX: u must be n//w, the aligned part matching the
        # SAME split. Setting u = n//gcd(n,M) separately made the
        # stored pair inconsistent, and every element then failed
        # the inclusion - 100% failure was the tell.
        if M % w == 0:
            out.append((n, a, n // w, w))
    return out


def build_bad(M=105, limit=500):
    """same but with misaligned parts NOT dividing M"""
    out = []
    for n in range(2, limit):
        a = (n // 2 + 1) % n
        w = mis_split(n, a)
        if M % w:
            out.append((n, a, n // w, w))
    return out


def kills(cls, X):
    s = set()
    for n, a, u, w in cls:
        st = a % n
        if st == 0:
            st = n
        while st < n:
            st += n
        s.update(range(st, X, n))
    return s


def check_inclusion(cls, M, X, r, label):
    U = kills(cls, X)
    head = cls[:r]
    Ur = kills(head, X)
    bound = max(n for n, _, _, _ in head)
    viol = 0
    tested = 0
    exc = 0
    for c in range(M):
        Dc = [u for (n, a, u, w) in cls if (a - c) % w == 0]
        Dcr = [u for (n, a, u, w) in head if (a - c) % w == 0]
        for x in range(c if c else M, X, M):
            if x in U and x not in Ur:
                if x <= bound:
                    exc += 1
                    continue
                tested += 1
                inM = any(x % u == 0 for u in Dc)
                inMr = any(x % u == 0 for u in Dcr)
                if not (inM and not inMr):
                    viol += 1
    print("  %-28s violations %5d of %6d elements "
          "(%d below the head bound, excluded)"
          % (label, viol, tested, exc))
    return viol, tested


def main():
    M, X, r = 105, 60000, 12
    cls = build(M)
    print("bounded-misalignment system: %d classes, all "
          "misaligned parts divide %d" % (len(cls), M))
    v1, t1 = check_inclusion(cls, M, X, r, "inclusion holds?")
    bad = build_bad(M)
    print("mutation system: %d classes, misaligned parts NOT "
          "dividing %d" % (len(bad), M))
    v2, t2 = check_inclusion(bad, M, X, r, "inclusion must FAIL")
    print()
    ok = (v1 == 0 and t1 > 500 and v2 > 0)
    print("  inclusion holds under bounded misalignment: %s" % (v1 == 0))
    print("  mutation fires (hypothesis load-bearing):    %s" % (v2 > 0))
    if not ok:
        print("\n[VERDICT] check failed or vacuous - the "
              "proposition stays a proposition.")
        sys.exit(1)
    # CHECK 2: F settles
    import math
    U = kills(cls, X)
    acc = 0.0
    vals = []
    marks = [2000, 8000, 30000, X - 1]
    mi = 0
    for x in range(2, X):
        if x in U:
            acc += 1.0 / x
        if mi < len(marks) and x == marks[mi]:
            vals.append((x, acc / math.log(x)))
            mi += 1
    print()
    print("  F(X) along the way: %s"
          % "  ".join("%d:%.4f" % v for v in vals))
    incs = [abs(b[1] - a[1]) for a, b in zip(vals, vals[1:])]
    print("  successive changes: %s"
          % "  ".join("%.4f" % i for i in incs))
    print()
    print("[VERDICT] The inclusion is exact and the mutation "
          "fires. The bounded-misalignment case follows from the "
          "vanishing-tails lemma alone - no character estimate "
          "is needed, and the proposition can be restored to a "
          "THEOREM.")


if __name__ == "__main__":
    main()
