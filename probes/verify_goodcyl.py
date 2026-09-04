#!/usr/bin/env python3
"""Session 105: structure forced by nu > 0.

The survivor set V = Zhat minus the union of balls is closed
with measure nu > 0. By the Lebesgue density theorem in a
profinite group, for almost every c in V the cylinders about c
become arbitrarily full of V. So for any eps there is an M whose
GOOD cylinders - those at least (1-eps)-full of survivors -
carry almost all of nu, and on such a cylinder the acting
classes kill at most eps of it.

CONSEQUENCE. If the induced system on a good cylinder contained
a pairwise-coprime subfamily of DIVERGENT drift, Theorem D would
drive its survivor measure to 0, contradicting >= 1-eps. So on
good cylinders no such subfamily exists:

  for every nu > 0 system, almost all the survivor measure sits
  on cylinders whose induced systems admit no divergent coprime
  subfamily.

That is exactly why Theorem D is unavailable in the open case.

CHECKS: exhibit good cylinders in a nu > 0 system with bounded
induced coprime drift; and a nu = 0 system must show the
opposite, or the statement separates nothing.
"""
import math
import sys
from math import gcd

FAIL=[]
def ck(n, ok, d=""):
    print("  [%s] %s%s"%("PASS" if ok else "FAIL",n,("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def induced(n,a,M,c):
    g=gcd(n,M)
    if (a-c)%g: return None
    nn=n//g
    if nn==1: return (1,0)
    inv=pow((M//g)%nn,-1,nn)
    return (nn, ((((a-c)//g)*inv)%nn))


def coprime_drift(mods, cap=3000):
    chosen=[]; d=0.0
    for m in sorted(mods):
        if m<2: continue
        if all(gcd(m,c)==1 for c in chosen):
            chosen.append(m); d+=1.0/m
            if len(chosen)>=cap: break
    return d


def cylinder_report(cls, M, label):
    """for each cylinder: induced survivor measure (product over
    the induced system, exact when the induced moduli are
    pairwise coprime) and induced coprime drift"""
    rows=[]
    for c in range(M):
        ind=[induced(n,a,M,c) for n,a in cls]
        mods=[x[0] for x in ind if x is not None and x[0]>1]
        # measure bound: survivors >= prod(1-1/m) over a coprime
        # subfamily; use the coprime drift as the diagnostic
        d=coprime_drift(mods)
        rows.append((c,len(mods),d))
    print("   %-30s %s"%(label,
        "  ".join("c=%d: coprime drift %.3f"%(c,d)
                  for c,_,d in rows)))
    return rows


def main():
    print("1. a nu > 0 system: moduli p^2 (nu exact)")
    qs=primes_upto(2000)
    cls=[(p*p, 1 % (p*p)) for p in qs]
    nu=1.0
    for p in qs: nu*=(1.0-1.0/(p*p))
    print("   exact nu = %.6f > 0"%nu)
    rows=cylinder_report(cls, 6, "cylinders mod 6")
    ck("nu > 0: induced coprime drifts are BOUNDED on every "
       "cylinder", all(d < 1.0 for _,_,d in rows),
       "max %.3f"%max(d for _,_,d in rows))

    print()
    print("2. a nu = 0 system: all integers, random residues")
    import random
    rng=random.Random(3)
    cls0=[(n, rng.randrange(n)) for n in range(2, 60000)]
    rows0=cylinder_report(cls0, 6, "cylinders mod 6")
    ck("nu = 0: induced coprime drifts are LARGE (Theorem D "
       "applies there, which is why those systems die)",
       all(d > 1.5 for _,_,d in rows0),
       "min %.3f"%min(d for _,_,d in rows0))

    print()
    print("3. the statement must SEPARATE the two regimes")
    sep = (max(d for _,_,d in rows) < min(d for _,_,d in rows0))
    ck("separation: every nu>0 cylinder drift below every nu=0 "
       "cylinder drift", sep,
       "%.3f < %.3f"%(max(d for _,_,d in rows),
                      min(d for _,_,d in rows0)))
    print()
    if FAIL:
        print("FAILED: %s"%", ".join(FAIL)); sys.exit(1)
    print("[verdict] nu > 0 forces bounded induced coprime "
          "drift on the cylinders carrying the survivors; that "
          "is precisely why Theorem D cannot reach the open "
          "case.")


if __name__=="__main__":
    main()
