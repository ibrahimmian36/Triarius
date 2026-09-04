#!/usr/bin/env python3
"""Session 100: measure the excursion set where nu is EXACT.

Every earlier attempt estimated nu from the same data as w,
which makes the comparison circular. Pairwise-coprime systems
break the circle: by CRT independence their survivor measure is
exactly prod(1 - 1/n_i), computed without touching the window.

System: moduli p^2 over primes (pairwise coprime, convergent
drift), where Theorem B already guarantees the logarithmic
density exists - so the excursion set MUST come out with
vanishing share, validating the measurement before it is
trusted.

CONTROL: a system whose excursion set is non-empty, or the
measurement cannot detect anything at all.

Not measured, and deliberately: the doubled construction kills
only even numbers, so w <= 1/2 = 1-nu identically and it can
never show an excursion. Knowing that in advance saves the run.
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


def wprofile(cls, X):
    m = bytearray(X)
    for n,a in cls:
        s=a%n
        if s==0: s=n
        while s<n: s+=n
        if s<X: m[s::n]=b"\x01"*len(m[s::n])
    out=[]; Y=4096
    while Y<=X:
        out.append((Y, sum(m[Y//2:Y])/float(Y-Y//2)))
        Y*=2
    return out


def main():
    X = 4*10**6
    print("1. coprime system with EXACT nu (Theorem B applies)")
    ps = primes_upto(2000)
    cls = [(p*p, 1 % (p*p)) for p in ps]
    nu = 1.0
    for p in ps:
        nu *= (1.0 - 1.0/(p*p))
    drift = sum(1.0/(p*p) for p in ps)
    print("   %d classes, drift %.5f (convergent), EXACT "
          "nu = %.6f, so 1-nu = %.6f"
          % (len(cls), drift, nu, 1-nu))
    w = wprofile(cls, X)
    print("   w(X): %s" % "  ".join("%.5f" % v for _,v in w))
    for delta in (0.01, 0.05):
        hits=[Y for Y,v in w if v > 1-nu+delta]
        print("   excursions with delta=%.2f: %d of %d windows"
              % (delta, len(hits), len(w)))
    ck("excursion set vanishes on a system Theorem B covers",
       all(v <= 1-nu+0.01 for _,v in w),
       "max w = %.5f vs 1-nu = %.5f"
       % (max(v for _,v in w), 1-nu))

    print()
    print("2. CONTROL: a system that MUST show excursions")
    # spike classes: moduli in (Y/2,Y] each killing one element
    spike=[]
    for Y in (262144, 1048576, 4194304):
        lo,hi=Y//2,Y; tgt=lo
        for n in range(lo+1,hi):
            if tgt>=hi: break
            spike.append((n,tgt%n)); tgt+=1
    spike=sorted(set(spike))
    # nu for THIS system is not a product; use the trivial
    # bound 1-nu <= drift, so any w above the drift is an
    # excursion beyond anything the measure can explain
    d2=sum(1.0/n for n,_ in spike)
    w2=wprofile(spike,X)
    print("   %d classes, drift %.4f" % (len(spike), d2))
    print("   w(X): %s" % "  ".join("%.4f" % v for _,v in w2))
    big=[v for _,v in w2 if v>0.3]
    ck("control shows large window densities (measurement can "
       "detect excursions)", len(big)>0,
       "%d windows above 0.3, max %.4f"
       % (len(big), max(v for _,v in w2)))
    print()
    if FAIL:
        print("FAILED: %s"%", ".join(FAIL)); sys.exit(1)
    print("[verdict] on a system where nu is exact and the "
          "theorem is known, the excursion set is empty; the "
          "control confirms the measurement can see excursions "
          "when they exist.")


if __name__=="__main__":
    main()
