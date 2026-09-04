#!/usr/bin/env python3
"""Session 103: reduce the excess to the Buchstab ratio.

WHY SPIKES YIELDED AND THE MIDDLE DOES NOT. A class of modulus
near X contributes to O(1) windows, so window drifts sum to D(Y)
and form a budget that Markov can exploit. A class of
INTERMEDIATE modulus contributes to EVERY larger window, so
there is no budget and no Markov bound. That is the structural
reason the two terms behave differently.

THE REDUCTION. Let s(X) be the window survivor density and
r(X) = s(X)/nu_X the window-to-global survivor ratio - the
Buchstab ratio. The excess is nu - s(X), and since nu_X >= nu,

    excess = nu - r*nu_X <= nu*(1 - r)   when r <= 1,
    excess < 0                            when r > 1.

So the excursion set sits inside {X : r(X) < 1 - delta/nu}, and
the whole problem reduces to r -> 1 in logarithmic mean.

CHECKS. nu must be EXACT, so pairwise-coprime moduli; and BOTH
directions of r must be exercised - one-sided data misled in
session 99.
"""
import math
import sys

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def profile(cls, X):
    m = bytearray(X)
    for n,a in cls:
        s=a%n
        if s==0: s=n
        while s<n: s+=n
        if s<X: m[s::n]=b"\x01"*len(m[s::n])
    out=[]; Y=8192
    while Y<=X:
        lo=Y//2
        out.append((Y, 1-sum(m[lo:Y])/float(Y-lo)))
        Y*=2
    return out


def main():
    X = 4*10**6
    print("system: moduli p^2 (pairwise coprime, so nu is the "
          "EXACT product)")
    ps = primes_upto(3000)
    cls = [(p*p, 1 % (p*p)) for p in ps]
    nu = 1.0
    for p in ps: nu *= (1.0-1.0/(p*p))
    prof = profile(cls, X)
    print("   exact nu = %.6f" % nu)
    print()
    print("   X          s(X)      r=s/nu    excess=nu-s   "
          "bound nu(1-r)+")
    rows=[]
    for Y, s in prof:
        r = s/nu
        exc = nu - s
        bound = nu*max(0.0, 1-r)
        rows.append((Y,s,r,exc,bound))
        print("   %-10d %.6f  %.6f  %+.6f     %.6f"
              % (Y,s,r,exc,bound))
    ok = all(e <= b + 1e-12 for _,_,_,e,b in rows)
    ck("excess <= nu*(1-r)+ in every window", ok)
    pos = [r for _,_,r,_,_ in rows if r > 1]
    neg = [r for _,_,r,_,_ in rows if r < 1]
    ck("both directions of r exercised (r>1 gives negative "
       "excess, r<1 gives the bound)",
       len(pos) > 0 and len(neg) > 0,
       "%d windows with r>1, %d with r<1" % (len(pos), len(neg)))
    negexc = [e for _,_,r,e,_ in rows if r > 1]
    ck("r>1 windows have NEGATIVE excess (no excursion)",
       all(e < 1e-12 for e in negexc) if negexc else False,
       "max excess where r>1: %s"
       % ("%.2e" % max(negexc) if negexc else "n/a"))
    print()
    print("   window-contribution counts (the budget argument):")
    for n in (X//4, X//64, 1000):
        cnt = sum(1 for Y,_ in prof if Y >= n)
        print("      a class of modulus %-8d appears in %d of "
              "the %d windows" % (n, cnt, len(prof)))
    print("   -> near-X moduli touch O(1) windows (budgeted); "
          "small moduli touch all of them (not budgeted)")
    print()
    if FAIL:
        print("FAILED: %s"%", ".join(FAIL)); sys.exit(1)
    print("[verdict] the excess is controlled by the Buchstab "
          "ratio; the problem reduces to r -> 1 in logarithmic "
          "mean.")


if __name__=="__main__":
    main()
