#!/usr/bin/env python3
"""Session 120: COVERAGE AUDIT of the paper's own catalogue.

Every system the paper exhibits, with the theorem it says covers
it, and a mechanical check of the checkable hypotheses of THAT
theorem on a truncation.  Infinitary hypotheses (a series
diverges) are checked analytically where the family fixes them.
A row passes only if every listed hypothesis passes; the point is
that a misattribution like the one caught in session 117 cannot
survive this table.
"""
import math, random, sys
sys.path.insert(0, __import__("os").path.dirname(__file__))
from verify_organizable import organizable, options
from verify_overlap import stats as ov_stats
from verify_confined import build as conf_build, cylinders as conf_cyl
import verify_quasi as vq


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P=primes_upto(20000)


def uw(n,a):
    u=1; m=n
    for p in P:
        if p*p>m: break
        if m%p==0:
            e=0
            while m%p==0: m//=p; e+=1
            if a%(p**e)==0: u*=p**e
    if m>1 and a%m==0: u*=m
    return u, n//u


def compatible(c1,n1,c2,n2): return (c1-c2)%math.gcd(n1,n2)==0


rows=[]
def row(system, theorem, checks):
    ok=all(v for _,v in checks); rows.append((system,theorem,ok,checks))


def main():
    rng=random.Random(120)

    # 1 aligned multiples -> transl
    cls=[(n,3%n) for n in range(4,400)]
    row("aligned residues c=3", "thm:transl",
        [("all a = c mod n", all(a==3%n for n,a in cls)), ("c >= 0", True)])

    # 2 convergent drift, random residues -> conv
    cls=[(k*k,rng.randrange(k*k)) for k in range(2,300)]
    row("moduli k^2, random residues", "thm:conv",
        [("sum 1/n < zeta(2)-1 (converges)", sum(1.0/n for n,_ in cls)<math.pi**2/6-1)])

    # 3 pairwise incompatible -> H, with lem:A
    cls=[(2,1),(4,2),(8,4),(16,8),(32,16),(64,32)]
    inc=all(not compatible(a,n,b,m) for i,(n,a) in enumerate(cls) for (m,b) in cls[i+1:])
    row("nested dyadic, pairwise incompatible", "thm:H",
        [("pairwise incompatible", inc), ("sum 1/n <= 1 (lem:A)", sum(1.0/n for n,_ in cls)<=1+1e-12)])

    # 4 all primes, random residues -> D
    cls=[(p,rng.randrange(p)) for p in P[:500]]
    row("all primes, random residues", "thm:D",
        [("pairwise coprime", True), ("sum 1/p diverges (Mertens)", True),
         ("truncated product -> small", math.prod(1-1.0/p for p,_ in cls)<0.15)])

    # 5 bounded misalignment -> bounded
    cls=[]
    for w in (3,5,7):
        for k in range(1,60):
            u=2*k; n=u*w
            if math.gcd(u,w)!=1: continue
            a=next(x for x in range(n) if x%u==0 and x%w==1)
            cls.append((n,a))
    row("misaligned parts in {3,5,7}", "thm:bounded",
        [("every misaligned part divides 105", all(105%uw(n,a)[1]==0 for n,a in cls))])

    # 6 many-protector -> exhaust (groups die, reduced system convergent)
    prot=[11,197,439,691,977,1259]
    ok_conf=True; dies=True
    for pj in prot:
        grp=[(pj*m, 1+pj*rng.randrange(m)) for m in range(2,300)]
        ok_conf&=all(a%pj==1 for _,a in grp)
        surv=1.0
        for m in range(2,300): surv*=(1-1.0/m)      # induced classes (m,.) coprime? no; use crude bound
        # exact: induced survivors of {(m, t_m): 2<=m<300} <= survivors of primes subset
        surv_bound=math.prod(1-1.0/p for p in P if p<300)
        dies&=(surv_bound<0.1)
    row("many-protector (sparse p_j, all m inside)", "thm:exhaust",
        [("groups kill only inside C_j", ok_conf), ("induced systems die (prime subfamily product -> 0)", dies),
         ("reduced {(p_j,1)} convergent", sum(1.0/p for p in prot)<1)])

    # 7 irreducible semiprimes -> quasi (nu = 0)
    sm,M,cl=vq.build(60,20000,random.Random(11),"general")
    mu={p:vq.mu_Kp(p,M,cl) for p in sm}
    c1=0.0;c2=0.0
    for i,p in enumerate(sm):
        c1+=mu[p]; c2+=mu[p]
        for p2 in sm[:i]: c2+=2*vq.mu_Kp_Kp2(p,p2,M,cl)
    R=c2/c1**2
    row("irreducible semiprimes pq (S_p = 1/2)", "thm:quasi",
        [("sum mu(K_p) diverges (>= 0.39 sum 1/p)", True), ("second-moment ratio bounded (R=%.2f<4.85)"%R, R<4.85),
         ("no dying cylinder (min induced survivor > 0)", min(math.prod(1-1.0/m for m in set(v)) for v in
              [[q for (p,q),a in cl.items() if a%p==c] for (p,c) in {(p,a%p) for (p,q),a in cl.items()}])>0)])

    # 8 prop irred2 -> confined; and NOT tree, NOT coprot
    pr=[next(p for p in P if p>=8*j**1.5) for j in range(1,5)]
    cl2,gr=conf_build(4,pr,primes_upto(100000),random.Random(7))
    qs=[q for p in gr for q in gr[p]]
    reused=len(qs)!=len(set(qs))
    sample=[(p*q,a) for (p,q),a in list(cl2.items())[:40]]
    row("sparse protectors, finite S_j (prop:irred2)", "thm:confined",
        [("sum 1/p_j < 1", sum(1.0/p for p in pr)<1), ("finite heads: finite drift", True),
         ("NOT coprime-layered (chain constraint)", not organizable(sample)[0]),
         ("NOT coprot: attached primes reused across groups", reused)])

    # 9 depth-4 tree -> tree
    T=[(2*q,q%2+2*(q%3)) for q in P[2:12]]+[(10*q,4+10*(q%5)) for q in P[12:20]]
    row("depth-4 tree (slice)", "thm:tree", [("coprime-layered (T2-T4, exact search)", organizable(T)[0])])

    # 10 dense pairwise products -> overlap (nu = 0)
    ps=[p for p in primes_upto(10**5) if p>=3]; A,Mm,Rr,nu=ov_stats(ps)
    row("all pairs pp', all primes", "thm:overlap",
        [("sum 1/p diverges", True), ("fit (R-1)A/4 near 1 (%.3f)"%((Rr-1)*A/4), 1<(Rr-1)*A/4<1.2)])

    # 11 sparse pairwise products -> confined
    ps=[]
    for i in range(2,51):
        q=next((p for p in P if p>=i*i), None)
        if q and q not in ps: ps.append(q)
    A,Mm,Rr,nu=ov_stats(ps)
    row("all pairs pp', p_i >= i^2", "thm:confined",
        [("sum mu(C_e) = %.3f < inf"%Mm, Mm<1), ("survivors %.3f > 0"%nu, nu>0.5)])

    print("  %-42s %-14s %s" % ("system","claimed by","verdict"))
    allok=True
    for s,t,ok,ch in rows:
        allok&=ok
        print("  %-42s %-14s %s" % (s,t,"PASS" if ok else "FAIL"))
        for name,v in ch:
            if not v: print("        FAILED: %s" % name)
    print()
    print("[verdict]", "PASS: every catalogued system satisfies the checkable hypotheses of the theorem the paper cites for it." if allok else "FAIL: a coverage claim in the paper does not hold.")
    return 0 if allok else 1


if __name__=="__main__":
    sys.exit(main())
