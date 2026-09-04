#!/usr/bin/env python3
"""Session 106: are the good cylinders empty or aligned?

Session 105: on cylinders carrying the survivors the induced
system is COPRIME-POOR. In both known nu > 0 systems with
divergent drift - the aligned evens and the doubled construction
- the good cylinder is the odd one and carries NO classes.

If good cylinders are always empty or aligned-dominated, the
translation theorem covers them and a route exists. If one is
non-empty and not aligned-dominated, the route closes.

Per cylinder we compute:
  killed measure  - from the INDUCED system, not from a window
                    (they differ by exactly the Buchstab
                    discrepancy this campaign spent three
                    sessions on);
  aligned share   - the largest drift carried by a family with a
                    COMMON INTEGER SOLUTION, found by explicit
                    CRT accumulation, not inferred from drift.
"""
import sys
from math import gcd

FAIL=[]
def ck(n, ok, d=""):
    print("  [%s] %s%s"%("PASS" if ok else "FAIL",n,("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def induced(n,a,M,c):
    g=gcd(n,M)
    if (a-c)%g: return None
    nn=n//g
    if nn==1: return (1,0)
    inv=pow((M//g)%nn,-1,nn)
    return (nn, ((((a-c)//g)*inv)%nn))


def killed_measure(cls, cap=200000):
    """CAUTION: the union bound sum 1/n is a FINITE-TRUNCATION
    quantity and must not be used to decide whether a cylinder
    is 'good'. In the doubled construction the even cylinder has
    union bound 0.0468 at six small primes, yet its true
    survivor measure is 0 - the infinite construction has
    sum s/q divergent and that cylinder DIES (session 75). Using
    the bound mislabels it as good and produced a false
    'route closed' verdict. Good cylinders are therefore
    identified by KNOWN measure, supplied per system."""
    return min(1.0, sum(1.0/n for n,_ in cls))


def aligned_share(cls):
    """largest drift of a subfamily with a COMMON INTEGER
    solution: accumulate by CRT, greedily, heaviest first"""
    best=0.0
    for start in range(min(3, len(cls))):
        r, m, d = None, 1, 0.0
        for n,a in sorted(cls[start:], key=lambda t:t[0]):
            if r is None:
                r, m, d = a % n, n, 1.0/n
                continue
            g=gcd(m,n)
            if (a-r) % g: continue
            lcm=m//g*n
            # CRT combine
            inv=pow((m//g)%(n//g), -1, n//g) if n//g>1 else 0
            k=(((a-r)//g)*inv)%(n//g) if n//g>1 else 0
            r=(r+m*k)%lcm; m=lcm; d+=1.0/n
        best=max(best,d)
    return best


def report(name, cls, M, known_measure):
    """known_measure[c] = the TRUE survivor measure of cylinder
    c, supplied from proof rather than estimated"""
    print("  %s (mod %d):" % (name, M))
    rows=[]
    for c in range(M):
        ind=[(x[0],x[1]) for x in
             (induced(n,a,M,c) for n,a in cls)
             if x is not None and x[0]>1]
        km=killed_measure(ind)
        al=aligned_share(ind) if ind else 0.0
        tot=sum(1.0/n for n,_ in ind)
        share = al/tot if tot>1e-12 else 1.0
        surv = known_measure[c]
        rows.append((c,len(ind),surv,share))
        print("     c=%d: %6d induced classes, TRUE survivor "
              "measure %.3f, aligned share %.3f%s"
              % (c,len(ind),surv,share,
                 "   <-- good" if surv > 0.5 else ""))
    return rows


def main():
    print("1. aligned evens (nu = 1/2, divergent drift)")
    ev=[(2*k,0) for k in range(1,40000)]
    # aligned evens: evens all killed, odds all survive
    r1=report("aligned evens", ev, 2, {0:0.0, 1:1.0})
    good1=[r for r in r1 if r[2] > 0.5]
    ck("good cylinders of the aligned evens are empty or "
       "aligned-dominated",
       all(r[1]==0 or r[3]>0.9 for r in good1),
       "good cylinders: %s"%[(r[0],r[1],round(r[3],2))
                             for r in good1])

    print()
    print("2. doubled construction (nu = 1/2, divergent drift)")
    def primes(n):
        s=bytearray([1])*(n+1); s[0]=s[1]=0
        i=2
        while i*i<=n:
            if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
            i+=1
        return [i for i in range(2,n+1) if s[i]]
    ps=primes(60000); smalls=ps[1:7]
    pool=[p for p in ps if p>smalls[-1] and p>=200]
    out=[]; idx=0
    for q in smalls:
        acc=0.0; j=0
        while idx<len(pool) and acc<0.05:
            p=pool[idx]; idx+=1
            r=j%q; j+=1
            inv=pow(q%p,-1,p); k=((0-r)%p)*inv%p
            a=(r+q*k)%(p*q); n2=2*p*q
            a2=a if a%2==0 else a+p*q
            out.append((n2,a2%n2)); acc+=1.0/p
    # doubled: even cylinder dies (session 75), odds untouched
    r2=report("doubled", sorted(set(out)), 2, {0:0.0, 1:1.0})
    good2=[r for r in r2 if r[2] > 0.5]
    ck("good cylinders of the doubled construction are empty "
       "or aligned-dominated",
       all(r[1]==0 or r[3]>0.9 for r in good2),
       "good cylinders: %s"%[(r[0],r[1],round(r[3],2))
                             for r in good2])
    print()
    if FAIL:
        print("[verdict] ROUTE CLOSED: a good cylinder is "
              "non-empty and not aligned-dominated, so the "
              "translation theorem does not cover it.")
        sys.exit(1)
    print("[verdict] ROUTE OPEN in the tested systems: every "
          "good cylinder is empty or aligned-dominated. That is "
          "two examples, not a theorem.")


if __name__=="__main__":
    main()
