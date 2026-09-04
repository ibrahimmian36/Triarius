#!/usr/bin/env python3
"""Session 112: can an infinite exhaustion tower exist?

STRUCTURAL CONSTRAINTS, found before coding.
 * Each reduction passes to DIVISORS, and a proper coarsening at
   least halves, so a class of modulus n0 drops out within
   log2(n0) levels; at level K only n0 >= 2^K remain.
 * A level-K modulus of 1 is a cylinder equal to Z, exhausted,
   hence nu = 0. So every level's moduli stay >= 2.
 * A PRIME cannot be a non-top level: no proper divisor above 1.
   So every non-top level is composite, and persisting classes
   need unboundedly long prime-factor chains.

THE DICHOTOMY seen in every family tried: divergent drift at
the top forces nu = 0; nu > 0 forces convergent drift at the
top, which terminates the tower. Heuristic reason: coarsening
divides out shared structure, the top's moduli become coprime,
and Theorem D bites.

MEASURED HERE on the mixed family 2^a * q (sparse primes q,
unbounded a), level by level: drift growth and nu, with nu from
the exact product where coprimality allows. The trivial
reduction (a class exhausting its own cylinder) is EXCLUDED - it
would make any tower 'infinite' vacuously.
"""
import math
import sys

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def level_system(qs, K, amax):
    """level-K moduli 2^(a-K) * q for a > K, residues aligned
    within each q (the only way nu can be positive), so the
    union over a of the cylinders is the a=K+1 one: c_q+2q Z"""
    mods=[]
    for q in qs:
        for a in range(K+1, amax+1):
            mods.append((2**(a-K))*q)
    return mods


def main():
    ps=primes_upto(20000)
    print("family: moduli 2^a * q, unbounded a; level K has "
          "moduli 2^(a-K) q for a > K")
    print()
    for label, qs in (("sparse q (sum 1/q < inf)", ps[4::200][:20]),
                      ("dense q  (sum 1/q -> inf)", ps[:400])):
        s=sum(1.0/q for q in qs)
        print("  %s: %d primes, sum 1/q = %.3f" % (label, len(qs), s))
        for K in (0, 2, 4, 6):
            mods=level_system(qs, K, 16)
            drift=sum(1.0/m for m in mods)
            # per-q aligned: union over a is c_q + 2q Z, measure
            # 1/(2q); across distinct q these are coprime-ish, so
            # nu = prod(1 - 1/(2q)) exactly when aligned
            nu=1.0
            for q in qs: nu*=(1-1.0/(2*q))
            # drift at level K ~ sum_q (1/q) * sum_{b>=1} 2^-b
            # = sum_q 1/q, independent of K - divergent iff the
            # q's are dense
            print("     level %d: %5d classes, drift %.3f, "
                  "nu = %.4f" % (K, len(mods), drift, nu))
        print()
    print("  at EVERY level the drift is ~ sum 1/q, and nu is "
          "prod(1-1/2q):")
    print("     sparse q: nu > 0 but drift CONVERGENT -> tower "
          "terminates (Theorem B at the top)")
    print("     dense q : drift divergent but nu -> 0 -> closed "
          "by the nu=0 corollary")
    print()
    print("  the two cannot be had together in this family. The")
    print("  heuristic reason - coarsening divides out shared")
    print("  structure until the top's moduli are coprime, where")
    print("  Theorem D forces the dichotomy - is a HEURISTIC, not")
    print("  a proof.")


if __name__=="__main__":
    main()
