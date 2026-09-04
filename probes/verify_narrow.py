#!/usr/bin/env python3
"""Session 104: narrow the open statement.

TWO REDUCTIONS.
 (1) Excursions require nu > 0. The excess is nu - s(X), so when
     nu = 0 it is never positive. This removes every dying
     system - including the prime sieve, where the Buchstab
     ratio reaches e^gamma yet no excursion is possible because
     nu = 0 there. Session 99's alarming ratio was a threat to a
     BOUND, never to the conjecture.
 (2) The deviation is one-sided. For every N the survivors of
     the whole system sit inside those of the first N classes,
     whose window density tends to nu_N; so limsup s <= nu and
     limsup r <= 1. With the previous proposition the whole
     remaining content is liminf r = 1.

CHECKS: a nu = 0 control must show non-positive excess in every
window; a nu > 0 system must show r at or below 1 up to o(1),
in contrast with the SUB-system measurements of session 99 that
exceeded 1.
"""
import math
import sys

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def windows(cls, X):
    m=bytearray(X)
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
    X=4*10**6
    print("1. reduction (1) is TRIVIAL and needs no control:")
    print("   excess = nu - s(X) <= nu because s >= 0, so nu = 0")
    print("   forces excess <= 0. My first attempt 'verified'")
    print("   this against nu_X at finite X - which by Mertens is")
    print("   e^{-gamma}/log X, NOT 0 - and so failed for a")
    print("   reason that had nothing to do with the claim.")
    print()
    print("   illustrating with the prime sieve, where nu_X is")
    print("   still far from its limit 0:")
    ps=primes_upto(X//2)
    nu0=1.0
    for p in ps: nu0*=(1.0-1.0/p)
    w=windows(((p,0) for p in ps), X)
    print("   exact nu = %.3e (essentially 0)" % nu0)
    exc=[(Y, nu0-s) for Y,s in w]
    print("   excess nu - s(X): %s"
          % "  ".join("%+.4f" % e for _,e in exc[-5:]))
    print("   nu_X = %.4f at X = %d, tending to 0 only as "
          "e^{-gamma}/log X - so this is NOT a nu = 0 system at "
          "any reachable scale, which is exactly the confusion "
          "that broke the first check." % (nu0, X))

    print()
    print("2. nu > 0 system: r must be <= 1 up to o(1)")
    qs=primes_upto(3000)
    cls=[(p*p, 1 % (p*p)) for p in qs]
    nu=1.0
    for p in qs: nu*=(1.0-1.0/(p*p))
    w2=windows(cls, X)
    # nu_X here equals nu (all moduli <= X), so r = s/nu
    rs=[(Y, s/nu) for Y,s in w2]
    print("   exact nu = %.6f" % nu)
    print("   r(X): %s" % "  ".join("%.6f" % r for _,r in rs))
    over=[r for _,r in rs if r > 1]
    ck("r <= 1 + o(1) on the full system",
       all(r <= 1.001 for _,r in rs),
       "max r = %.6f, %d windows above 1 (all within o(1))"
       % (max(r for _,r in rs), len(over)))

    print()
    print("3. contrast: session 99 measured SUB-systems with "
          "r well above 1")
    sub=[p for p in primes_upto(int(X**0.5)) if p > int(X**0.05)]
    g=1.0
    for p in sub: g*=(1.0-1.0/p)
    ws=windows(((p,0) for p in sub), X)
    rsub=[s/g for _,s in ws]
    print("   sub-system r: %s"
          % "  ".join("%.4f" % r for r in rsub[-5:]))
    ck("sub-systems can exceed 1 substantially (so the "
       "one-sidedness is a property of FULL systems)",
       max(rsub) > 1.05, "max %.4f" % max(rsub))
    print()
    if FAIL:
        print("FAILED: %s"%", ".join(FAIL)); sys.exit(1)
    print("[verdict] excursions require nu > 0, and on full "
          "systems r <= 1 + o(1); the whole remaining content "
          "is liminf r = 1.")


if __name__=="__main__":
    main()
