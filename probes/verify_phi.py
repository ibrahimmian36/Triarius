#!/usr/bin/env python3
"""Session 84 self-audit: Theorem 5.2's proof appeals to
Davenport-Erdos 'on each cylinder', but D-E is a theorem about
sets of MULTIPLES in N, and a set of multiples intersected with
a fixed residue class mod M is not a set of multiples: writing
x = c + My turns divisibility into a congruence, so the step is
unjustified as written.

The repair is the MONOTONE-WEIGHT generalisation of D-E's
Lemma 1. For a system whose misaligned parts divide M, define

   phi(x) = (1/M) * #{ c mod M :
                       x = a_i (mod w_i) for some class i
                       applicable at x },

where class i is APPLICABLE at x when u_i | x and x >= n_i.
Because the applicable set only GROWS when x is multiplied,
phi is monotone under divisibility, and Lemma 1's one-line proof
goes through verbatim for any [0,1]-valued divisibility-monotone
weight.

CHECK: phi(d) <= phi(n) whenever d | n, exactly, over a
bounded-misaligned system.
MUTATION: a weight built from a NON-divisibility-closed
applicability rule (u_i does NOT divide x) must violate
monotonicity - otherwise the check is vacuous.
"""
from math import gcd

def mis_split(n, a):
    """n = u*w, gcd(u,w)=1, u | a; smallest such w"""
    best = None
    for w in range(1, n + 1):
        if n % w:
            continue
        u = n // w
        if gcd(u, w) == 1 and a % u == 0:
            if best is None or w < best:
                best = w
    return best


def build(limit=60):
    """A system in which NO class is aligned. Any class with
    misaligned part w = 1 has a vacuous cylinder condition and
    single-handedly forces phi == 1, which is why the first two
    attempts saturated and their mutations could not fire. Here
    every misaligned part is 3, 5 or 7 - the same classes used to
    verify the cylinder reduction of Theorem 5.2."""
    base = [(6, 2), (10, 4), (14, 6), (15, 5), (21, 7)]
    out = []
    for n, a in base:
        w = mis_split(n, a)
        assert w > 1, (n, a, w)
        out.append((n, a, n // w, w))
    return out


def phi(x, cls, M, applicable):
    hit = set()
    for (n, a, u, w) in cls:
        if not applicable(x, n, u):
            continue
        for c in range(M):
            if (c - a) % w == 0:
                hit.add(c)
    return len(hit) / float(M)


def run(name, applicable, cls, M=105, X=1200):
    viol = 0
    tested = 0
    worst = None
    cache = {}
    for n in range(2, X):
        for d in range(2, n):
            if n % d:
                continue
            tested += 1
            if d not in cache:
                cache[d] = phi(d, cls, M, applicable)
            if n not in cache:
                cache[n] = phi(n, cls, M, applicable)
            if cache[d] > cache[n] + 1e-12:
                viol += 1
                if worst is None:
                    worst = (d, n, cache[d], cache[n])
    print("  %-34s violations %5d of %6d divisor pairs%s"
          % (name, viol, tested,
             ("   witness phi(%d)=%.4f > phi(%d)=%.4f" % worst)
             if worst else ""))
    return viol, tested


def main():
    cls = build()
    M = 105
    vals = sorted({round(phi(x, cls, M,
                             lambda x, n, u: x % u == 0
                             and x >= n), 4)
                   for x in range(2, 400)})
    print("bounded-misaligned system: %d classes, misaligned "
          "parts all divide 105" % len(cls))
    print("phi takes %d distinct values in [%.3f, %.3f] -- it "
          "must NOT be saturated at 1, or the test is vacuous"
          % (len(vals), vals[0], vals[-1]))
    if len(vals) < 3 or vals[-1] == vals[0]:
        print("VACUOUS: phi is constant; shrink the system further")
        return
    print()
    v1, t1 = run("phi, divisibility applicability",
                 lambda x, n, u: x % u == 0 and x >= n, cls)
    v2, t2 = run("MUTATION: non-divisibility rule",
                 lambda x, n, u: x % u != 0 and x >= n, cls)
    print()
    ok = (v1 == 0 and t1 > 1000 and v2 > 0)
    print("phi is divisibility-monotone: %s" % (v1 == 0))
    print("mutation fires (non-vacuous):  %s" % (v2 > 0))
    print("-> %s" % ("monotone-weight route is SOUND; "
                     "Theorem 5.2's proof should use it rather "
                     "than appeal to D-E cylinderwise"
                     if ok else "CHECK FAILED - weaken the "
                     "theorem rather than patch it"))


if __name__ == "__main__":
    main()
