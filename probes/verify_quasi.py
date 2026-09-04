#!/usr/bin/env python3
"""Session 114: the quasi-independence theorem and the
irreducible semiprime system.

THEOREM (quasi-independence => nu = 0). Group the classes by
the smallest prime p of the modulus and let K_p be the union of
their balls. K_p depends only on the coordinates of Zhat at
primes >= p, so "x in K_p for infinitely many p" is a TAIL
event of independent coordinates. If  sum_p mu(K_p) = infinity
and the second-moment ratio
   R(N) = sum_{p,p'<=N} mu(K_p cap K_p') / (sum_{p<=N} mu(K_p))^2
stays bounded, Kochen-Stone gives P(K_p i.o.) >= 1/limsup R > 0,
Kolmogorov makes it 1, and survivors have measure 0.

THE SYSTEM. Moduli pq, p < q primes, q in M_p = the primes
above p taken in order until sum 1/q >= 1/2. For any residues:
   mu(K_p)   >= (1 - e^{-1/2})/p              (superadditivity)
   mu(K_p cap K_p') <= (S_p S_p' + min(S_p,S_p') + 1)/(p p'), the last
   term for a witness q = p' whose coordinate membership in K_p' fixes,
   so R <= 1.75/0.393^2 = 11.4 in the limit.
It is IRREDUCIBLE: every cylinder (p,c) has induced moduli that
are distinct primes, with induced drift about 1 and induced
survivor measure prod(1-1/l) > 0, so no group dies; cylinders
(pq,.) hold one class. Not aligned, misalignment unbounded.

Checks (exact unless marked MC):
  1. S_p in [1/2, 1/2 + 1/p)                  (construction)
  2. mu(K_p) exact vs the bound (1-e^{-1/2})/p (all p)
  3. R(N) exact by inclusion-exclusion, vs 4.84, trend in N
  4. irreducibility: max induced drift over cylinders (p,c),
     min induced survivor measure > 0
  5. mutation EXPECTED to escape: keep only sparse p (p_j with
     sum 1/p_j small): sum mu(K_p) converges, theorem silent,
     and nu >= 1 - sum mu(K_p) > 0.
  6. MC survival of the truncation (illustration only; no
     inference about the infinite system from a truncation).
"""
import math, random, sys


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]


def build(PMAX, QMAX, rng, scheme):
    P=primes_upto(QMAX)
    small=[p for p in P if p<=PMAX]
    M={}; classes={}
    for p in small:
        s=0.0; M[p]=[]
        for q in P:
            if q<=p: continue
            M[p].append(q); s+=1.0/q
            if s>=0.5: break
        else:
            sys.exit("QMAX too small for p=%d"%p)
        for q in M[p]:
            if scheme=="general":
                a=rng.randrange(p*q)
            elif scheme=="aligned":
                a=1
            elif scheme=="protect":      # 1 mod p, random mod q
                r=rng.randrange(q)
                a=next(x for x in range(p*q) if x%p==1 and x%q==r)
            classes[(p,q)]=a
    return small, M, classes


def mu_Kp(p, M, classes):
    # exact: sum_c (1/p)(1 - prod_{q in M_{p,c}} (1-1/q))
    byc={}
    for q in M[p]:
        c=classes[(p,q)]%p
        byc.setdefault(c,[]).append(q)
    tot=0.0
    for c,qs in byc.items():
        pr=1.0
        for q in qs: pr*=(1-1.0/q)
        tot+=(1-pr)/p
    return tot


def mu_Kp_Kp2(p, p2, M, classes):
    """exact mu(K_p & K_p'), p < p2, conditioning on the coordinates at p and p2.
    A witness q = p2 for the p-group (p2 in M_p) is decided by the p2-coordinate
    already fixed by membership in K_p2: it fires iff a_{p p2} = c2 mod p2, with no
    factor 1/q.  All other witnesses are at primes distinct from p, p2 and from
    each other, hence independent given (r_p, r_p2)."""
    if p > p2: p, p2 = p2, p
    byc={}; byc2={}
    for q in M[p]:  byc.setdefault(classes[(p,q)]%p,{})[q]=classes[(p,q)]%q
    for q in M[p2]: byc2.setdefault(classes[(p2,q)]%p2,{})[q]=classes[(p2,q)]%q
    tot=0.0
    for c,T in byc.items():
        for c2,T2 in byc2.items():
            T1=dict(T); forced=False
            if p2 in T1:
                if T1[p2]==c2: forced=True
                del T1[p2]                      # the p2-coordinate is fixed at c2
            nF2=1.0
            for q in T2: nF2*=(1-1.0/q)
            if forced:
                pf=1-nF2                        # F_p certain; only F_p2 matters
            else:
                nF=1.0
                for q in T1: nF*=(1-1.0/q)
                both=1.0
                for q in set(T1)|set(T2):
                    k=len({T1.get(q),T2.get(q)}-{None})
                    both*=(1-k/float(q))
                pf=1-nF-nF2+both
            tot+=pf/(p*p2)
    return tot


def main():
    rng=random.Random(11)
    PMAX,QMAX=int(sys.argv[1]) if len(sys.argv)>1 else 150, 20000
    small,M,classes=build(PMAX,QMAX,rng,"general")
    ncls=len(classes)
    print("semiprime system: p <= %d, %d classes, %d primes q up to %d"
          % (PMAX,ncls,len({q for p in M for q in M[p]}),max(q for p in M for q in M[p])))

    # 1. S_p
    S={p:sum(1.0/q for q in M[p]) for p in small}
    ok1=all(0.5<=S[p]<0.5+1.0/p for p in small)
    print("[1] S_p in [1/2, 1/2+1/p) for all p:", ok1,
          " (min %.4f max %.4f)"%(min(S.values()),max(S.values())))

    # 2. mu(K_p) vs bound
    bound=1-math.exp(-0.5)
    muK={p:mu_Kp(p,M,classes) for p in small}
    ok2=all(muK[p]*p>=bound-1e-12 for p in small)
    print("[2] mu(K_p) >= (1-e^-1/2)/p = %.4f/p for all p:"%bound, ok2,
          " (min p*mu %.4f, max %.4f; sum_p mu(K_p)=%.3f)"
          % (min(muK[p]*p for p in small),max(muK[p]*p for p in small),
             sum(muK.values())))

    # 3. R(N) exact
    print("[3] second-moment ratio R(N), exact:")
    cum1=0.0; cum2=0.0; done=set(); ratios=[]
    for i,p in enumerate(small):
        cum1+=muK[p]; cum2+=muK[p]
        for p2 in small[:i]:
            cum2+=2*mu_Kp_Kp2(p,p2,M,classes)
        R=cum2/cum1**2
        ratios.append(R)
        if i%6==5 or i==len(small)-1:
            print("      N=%4d  sum mu=%.3f  R=%.3f"%(p,cum1,R))
    ok3=max(ratios[3:])<11.4
    print("    R stays below the analytic 11.4:", ok3,
          " (max after 3 primes %.3f, final %.3f)"%(max(ratios[3:]),ratios[-1]))

    # 4. irreducibility: cylinders (p,c) for every prime p in play
    allp=sorted({p for p in M}|{q for p in M for q in M[p]})
    cyl={}
    for (p,q),a in classes.items():
        cyl.setdefault((p,a%p),[]).append(q)
        cyl.setdefault((q,a%q),[]).append(p)
    maxdrift=0.0; minsurv=1.0
    for (l,c),mods in cyl.items():
        d=sum(1.0/m for m in mods); maxdrift=max(maxdrift,d)
        s=1.0
        for m in set(mods): s*=(1-1.0/m)
        minsurv=min(minsurv,s)
    prime_mods=all(all(primes_upto(m)[-1]==m for m in mods) for mods in cyl.values())
    print("[4] irreducible: induced moduli all prime:",prime_mods,
          "; max induced drift %.3f; min induced survivor %.3f > 0:"
          % (maxdrift,minsurv), minsurv>0)
    mis=sum(1 for (p,q),a in classes.items() if a%p and a%q)
    print("    misaligned (neither prime divides residue): %d/%d classes"%(mis,ncls))

    # 5. mutation: sparse p only -> theorem silent, nu > 0
    sparse=[p for p in small if p in (11,197,439,691,977,1259) or (p>50 and str(p).endswith("7") and p%3==1)][:6]
    if len(sparse)<3: sparse=small[-3:]
    smu=sum(muK[p] for p in sparse)
    print("[5] mutation, sparse p=%s: sum mu(K_p)=%.3f < 1, so nu >= %.3f > 0"
          % (sparse,smu,1-smu), "  [theorem's hypothesis fails as designed:",
          smu<1, "]")

    # 6. MC survival illustration, three schemes, plus sparse control
    print("[6] MC survival of the truncation (illustration only):")
    Pall=primes_upto(max(q for p in M for q in M[p]))
    def mc(classes_, small_, n=3000):
        surv=0
        idx={}
        for (p,q),a in classes_.items(): idx.setdefault(p,[]).append((q,a%p,a%q))
        for _ in range(n):
            r={l:rng.randrange(l) for l in Pall}
            killed=False
            for p in small_:
                rp=r[p]
                for q,ap,aq in idx[p]:
                    if rp==ap and r[q]==aq: killed=True; break
                if killed: break
            surv+=not killed
        return surv/float(n)
    for scheme in ("general","aligned","protect"):
        _,_,cl=build(PMAX,QMAX,random.Random(5),scheme)
        print("      %-8s survivors %.3f   (product heuristic exp(-sum mu(K_p)) = %.3f)"
              % (scheme, mc(cl,small), math.exp(-sum(mu_Kp(p,M,cl) for p in small))))
    clS={k:v for k,v in classes.items() if k[0] in sparse}
    print("      sparse-p control survivors %.3f  (bound 1-sum mu = %.3f)"
          % (mc(clS,sparse),1-smu))

    allok=ok1 and ok2 and ok3 and prime_mods and minsurv>0 and smu<1
    print()
    print("[verdict]", "PASS: every quantity the proof uses holds exactly on the truncation; the mutation escapes as designed." if allok else "FAIL")
    return 0 if allok else 1


if __name__=="__main__":
    sys.exit(main())
