#!/usr/bin/env python3
"""FINAL verification of every numerical claim in the paper.
Each claim is checked independently; several carry controls that
must fail if the property being checked were absent."""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd, log

FAILED = []
def ck(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                          ("  " + detail) if detail else ""))
    if not ok: FAILED.append(name)


def dens(cls, X):
    k = bytearray(X)
    for n, a in cls:
        s = a % n
        if s == 0: s = n
        while s < n: s += n
        k[s::n] = b"\x01" * len(k[s::n])
    lo = X // 2
    return 1.0 - sum(k[lo:X]) / float(X - lo)


print("1. Lemma A: exact maximum over incompatible families")
U=[(m,b) for m in range(2,15) for b in range(m)]
n=len(U)
inc=[[((U[i][1]-U[j][1])%gcd(U[i][0],U[j][0]))!=0 for j in range(n)] for i in range(n)]
order=sorted(range(n),key=lambda i:U[i][0])
w=[Fraction(1,U[i][0]) for i in order]
suf=[Fraction(0)]*(n+1)
for i in range(n-1,-1,-1): suf[i]=suf[i+1]+w[i]
best=[Fraction(0)]
def rec(i,ch,acc):
    if acc>best[0]: best[0]=acc
    if i>=n or acc+suf[i]<=best[0]: return
    v=order[i]
    if all(inc[v][u] for u in ch):
        ch.append(v); rec(i+1,ch,acc+w[i]); ch.pop()
    rec(i+1,ch,acc)
rec(0,[],Fraction(0))
ck("exact max = 1 over 104 classes (mod<=14)", best[0]==1 and n==104,
   "max=%s classes=%d"%(best[0],n))

fam=[(2,0)]+[(12,r) for r in (1,3,5,7,9,11)]
ck("extremal family incompatible, drift 1",
   sum(Fraction(1,m) for m,_ in fam)==1 and
   all((b1-b2)%gcd(m1,m2)!=0 for (m1,b1),(m2,b2) in combinations(fam,2)))

print("2. Theorems H and D")
X=600000
h=[(2**j,2**(j-1)) for j in range(1,15)]
ck("H: density = 1 - sum 1/n",
   abs(dens(h,X)-(1-sum(1.0/m for m,_ in h)))<0.002)
ps=[p for p in range(2,200) if all(p%q for q in range(2,int(p**.5)+1))]
cp=[(p,1%p) for p in ps]; pr=1.0
for m,_ in cp: pr*=(1-1.0/m)
ck("D: density = prod(1-1/n)", abs(dens(cp,X)-pr)<0.005)
mut=[(2*p,2%(2*p)) for p in ps[1:14]]; pm=1.0
for m,_ in mut: pm*=(1-1.0/m)
ck("  control: non-coprime family deviates", abs(dens(mut,X)-pm)>0.02,
   "deviation %.4f"%abs(dens(mut,X)-pm))

print("3. Proposition (the break), moduli 2..39, a=1")
killed=set()
for nn in range(2,40):
    s=1%nn
    if s==0: s=nn
    while s<nn: s+=nn
    killed.update(range(s,3000,nn))
ck("31 killed, 62 not", 31 in killed and 62 not in killed)
def Lam(k):
    for p in range(2,k+1):
        if k%p==0:
            m=k
            while m%p==0: m//=p
            return log(p) if m==1 else 0.0
    return 0.0
rhs=sum(Lam(42//d) for d in range(1,43) if 42%d==0 and d in killed)
ck("Lemma 1 fails at 42: lhs 0, rhs log2+log3+log7",
   42 not in killed and abs(rhs-(log(2)+log(3)+log(7)))<1e-9,
   "rhs=%.4f"%rhs)

print("4. Theorem (bounded misalignment): cylinder reduction")
def mis(nn,a):
    b=None
    for ww in range(1,nn+1):
        if nn%ww: continue
        uu=nn//ww
        if gcd(uu,ww)==1 and a%uu==0 and (b is None or ww<b): b=ww
    return b
cls=[(6,2),(10,4),(14,6),(15,5),(21,7)]
ws=[mis(nn,a) for nn,a in cls]; M=1
for ww in ws: M=M*ww//gcd(M,ww)
bad=tot=0
for c in range(M):
    for (nn,a),ww in zip(cls,ws):
        uu=nn//ww
        S={x for x in range(1,40000) if x%M==c and x%nn==a%nn}
        T={x for x in range(1,40000) if x%M==c and x%uu==0}
        if not S: continue
        tot+=1
        if S!=T: bad+=1
ck("cylinder reduction exact, w=3,5,7, M=105", bad==0 and tot>20,
   "%d/%d pairs exact"%(tot-bad,tot))

print("5. Remark (why phi does not bridge)")
info=[(nn,a,nn//mis(nn,a),mis(nn,a)) for nn,a in cls]
XX=200000
km=bytearray(XX+1)
for nn,a,uu,ww in info:
    s=a%nn
    if s==0: s=nn
    while s<nn: s+=nn
    for i in range(s,XX+1,nn): km[i]=1
F=sum(1.0/x for x in range(2,XX+1) if km[x])/log(XX)
def phi(x):
    hit=set()
    for nn,a,uu,ww in info:
        if x%uu==0 and x>=nn:
            for c in range(M):
                if (c-a)%ww==0: hit.add(c)
    return len(hit)/float(M)
# exact limits: killed set periodic mod 210; phi (no threshold) periodic mod 70
from fractions import Fraction as _Fr
_kill=sum(1 for x in range(210) if any(x%nn==a for nn,a,_,_ in info))
def _phi0(x):
    hit=set()
    for nn,a,uu,ww in info:
        if x%uu==0:
            for c in range(M):
                if (c-a)%ww==0: hit.add(c)
    return _Fr(len(hit),M)
_pm=sum(_phi0(x) for x in range(70))/70
ck("phi remark: lim F = 73/210 = 0.3476, log-mean of phi = 169/490 = 0.3449, different",
   _Fr(_kill,210)==_Fr(73,210) and _pm==_Fr(169,490) and _pm!=_Fr(73,210),
   "killed %s phi-mean %s"%(_Fr(_kill,210),_pm))

print("6. Refuted-routes section claims")
ck("HR extreme case: density 0.125 vs bound 0.328125",
   abs(dens([(2,0),(4,1),(8,3)],60000)-0.125)<0.002 and
   abs((1-.5)*(1-.25)*(1-.125)-0.328125)<1e-12)
ck("HR enumeration count = 53130",
   sum(a*b*c for a,b,c in combinations(range(2,13),3))==53130)
# block [N,2N] on Zhat: general position ~ product ~ 1/2, aligned falls
import random as _rb
def _blockmeasure(N, aligned, seed=5, T=20000):
    rng=_rb.Random(seed); P=[q for q in range(2,2*N+1) if all(q%d for d in range(2,int(q**.5)+1))]
    def fac(n):
        f={}; m=n
        for q in P:
            if q*q>m: break
            while m%q==0: f[q]=f.get(q,0)+1; m//=q
        if m>1: f[m]=f.get(m,0)+1
        return f
    mods=list(range(N,2*N+1)); F={n:fac(n) for n in mods}
    pe={q:max(f.get(q,0) for f in F.values()) for q in P}
    res={n:(0 if aligned else rng.randrange(n)) for n in mods}
    k=0
    for _ in range(T):
        r={q:rng.randrange(q**pe[q]) for q in P if pe[q]>0}
        k+=any(all(r[q]%(q**e)==res[n]%(q**e) for q,e in F[n].items()) for n in mods)
    return k/float(T)
_g40=_blockmeasure(40,False); _a40=_blockmeasure(40,True)
ck("block [40,80]: general position 0.513, aligned 0.413, aligned below general",
   round(_g40,3)==0.513 and round(_a40,3)==0.413 and _a40<_g40, "general %.3f aligned %.3f"%(_g40,_a40))
_gb={N:_blockmeasure(N,False) for N in (20,80)}; _ab={N:_blockmeasure(N,True) for N in (20,80)}
ck("block [N,2N] Haar measure on Zhat (MC over residues mod the relevant prime powers) at N=20,80: general 0.516, 0.507; aligned 0.426, 0.392",
   round(_gb[20],3)==0.516 and round(_gb[80],3)==0.507 and round(_ab[20],3)==0.426 and round(_ab[80],3)==0.392,
   "general %.3f %.3f aligned %.3f %.3f"%(_gb[20],_gb[80],_ab[20],_ab[80]))

# transient spikes: classes 0 mod n, n in (Y/2,Y], Y = 2^13, 2^16, 2^19, 2^21
import math as _sm
_Ys=[2**13,2**16,2**19,2**21]; _M=_Ys[-1]; _kl=bytearray(_M+1)
for _Y in _Ys:
    for _n in range(_Y//2+1,_Y+1): _kl[_n::_n]=b'\x01'*len(range(_n,_M+1,_n))
_wd=[];_F=[];_H=0.0
for _k in range(12,21):
    _lo,_hi=2**_k,2**(_k+1); _wd.append(sum(_kl[_lo+1:_hi+1])/(_hi-_lo))
    _H+=sum(1.0/_x for _x in range(_lo+1,_hi+1) if _kl[_x]); _F.append(_H/_sm.log(_hi))
ck("spike display: window killed 1.000 0.583 0.527 1.000 0.692 0.645 1.000 0.747 1.000",
   [round(v,3) for v in _wd]==[1.0,0.583,0.527,1.0,0.692,0.645,1.0,0.747,1.0], str([round(v,3) for v in _wd]))
ck("spike display: F rises 0.077 0.112 0.140 0.194 0.222 0.246 0.286 0.308 0.341, monotone",
   [round(v,3) for v in _F]==[0.077,0.112,0.140,0.194,0.222,0.246,0.286,0.308,0.341] and all(_F[i]<_F[i+1] for i in range(8)),
   str([round(v,3) for v in _F]))

print("7. Refuted routes: transient spikes do not move the log average")
import math as _m
def _sieve(cls,X):
    m=bytearray(X)
    for n,a in cls:
        s=a%n
        if s==0: s=n
        while s<n: s+=n
        if s<X: m[s::n]=b"\x01"*len(m[s::n])
    return m
X2=4*10**6
cls2=[]
for Y in (8192,65536,524288,2097152):
    lo,hi=Y//2,Y; tgt=lo
    for nn in range(lo+1,hi):
        if tgt<nn: tgt=nn
        if tgt>=hi: break
        cls2.append((nn,tgt%nn)); tgt+=1
cls2=sorted(set(cls2))
m2=_sieve(cls2,X2)
ws=[]; Y=2048
while Y<X2:
    ws.append(sum(m2[Y//2:Y])/float(Y-Y//2)); Y*=2
acc=0.0; Fs=[]; Y=2048; i=0
marks=[]
Y=2048
while Y<X2: marks.append(Y); Y*=2
for x in range(2,X2):
    if m2[x]: acc+=1.0/x
    if i<len(marks) and x==marks[i]: Fs.append(acc/_m.log(x)); i+=1
osc=max(ws)-min(ws[2:])
Fmono=all(b>=a-1e-9 for a,b in zip(Fs,Fs[1:]))
ck("window density oscillates by >0.4 while F is monotone",
   osc>0.4 and Fmono, "oscillation %.3f, F monotone %s"%(osc,Fmono))

# --- Proposition (irreducible semiprimes), p <= 150 ---
import os, random as _r
_here=os.path.dirname(os.path.abspath(__file__))
for _c in (os.path.join(_here,"..","probes","erdos25"),
           os.path.join(_here,"..","probes")):
    if os.path.exists(os.path.join(_c,"verify_quasi.py")):
        sys.path.insert(0,_c); break
import verify_quasi as vq
_sm,_M,_cl=vq.build(150,20000,_r.Random(11),"general")
ck("irreducible semiprimes: 7,174 classes on 582 primes up to 4,253",
   len(_cl)==7174 and len({q for p in _M for q in _M[p]})==582
   and max(q for p in _M for q in _M[p])==4253)
_S={p:sum(1.0/q for q in _M[p]) for p in _sm}
ck("S_p in [1/2,1/2+1/p)", all(0.5<=_S[p]<0.5+1.0/p for p in _sm))
_mu={p:vq.mu_Kp(p,_M,_cl) for p in _sm}
_minpmu=min(_mu[p]*p for p in _sm)
ck("min p*mu(K_p) = 0.4815 >= 0.3935", round(_minpmu,4)==0.4815 and _minpmu>=1-_m.exp(-0.5))
_c1=0.0;_c2=0.0;_R=[]
for _i,p in enumerate(_sm):
    _c1+=_mu[p];_c2+=_mu[p]
    for p2 in _sm[:_i]: _c2+=2*vq.mu_Kp_Kp2(p,p2,_M,_cl)
    _R.append(_c2/_c1**2)
ck("ratio 1.889 at N=149, falling from 2.159 at N=13, below 11.4",
   round(_R[-1],3)==1.889 and round(_R[5],3)==2.159 and max(_R[3:])<11.4,
   "final %.3f"%_R[-1])
_cyl={}
for (p,q),a in _cl.items():
    _cyl.setdefault((p,a%p),[]).append(q); _cyl.setdefault((q,a%q),[]).append(p)
_md=max(sum(1.0/m for m in v) for v in _cyl.values())
_ms=min(__import__("functools").reduce(lambda x,m:x*(1-1.0/m),set(v),1.0) for v in _cyl.values())
ck("max induced drift 0.886, min induced survivor 0.316 > 0",
   round(_md,3)==0.886 and round(_ms,3)==0.316 and _ms>0)
_sp=sum(_mu[p] for p in (11,67,97,127))
ck("sparse-p control: sum mu = 0.061, survivors >= 0.939",
   round(_sp,3)==0.061 and round(1-_sp,3)==0.939)

# --- Proposition (irreducible with survivors), six groups ---
import verify_confined as vc
_P=vc.primes_upto(100000)
import math as _mth
_prot=[_P[_mth.ceil(4*j**1.5)-1] for j in range(1,7)]
_cl,_gr=vc.build(6,_prot,_P,_r.Random(7))
ck("irreducible with survivors: protectors 7..277 (the ceil(4j^1.5)-th primes), 6,500 classes",
   _prot==[7,37,73,131,197,277] and len(_cl)==6500)
_sp=sum(1.0/p for p in _prot)
_tail=sum(1.0/(k*_mth.log(k)) for k in (_mth.ceil(4*j**1.5) for j in range(7,200000)))
ck("sum 1/p_j = 0.1999 over six terms, whole sum < 0.24 (tail via p_k > k log k)", round(_sp,4)==0.1999 and _sp+_tail<0.24)
_cy=vc.cylinders(_cl); _Ps=set(_P); _ms=1.0; _md=0.0; _ap=True
for (l,c),mods in _cy.items():
    if len(mods)<2: continue
    _ap&=(l in _Ps) and all(m in _Ps for m in mods) and len(set(mods))==len(mods)
    _md=max(_md,sum(1.0/m for m in mods))
    _ms=min(_ms,__import__("functools").reduce(lambda x,m:x*(1-1.0/m),mods,1.0))
ck("every multi-class cylinder prime with distinct prime induced moduli; min survivor 0.0637; max drift 2.449",
   _ap and round(_ms,4)==0.0637 and round(_md,3)==2.449)
ck("6,490 misaligned classes", sum(1 for (p,q),a in _cl.items() if a%p and a%q)==6490)

# --- Theorem (coprime-layered trees): explicit tree of depth four ---
import verify_tree as vt
_nodes,_cells,_ms=vt.build(4,_r.Random(19))
_cls=[cl for n in _nodes for cl in n[4]]
ck("tree: 617 nodes, 23,361 classes", len(_nodes)==617 and len(_cls)==23361)
_nu=0.0
for d,c,br in _cells:
    _s=1.0
    for q in br: _s*=(1-1.0/q)
    _nu+=_s/d
ck("tree survivor measure 0.96432 by branch product", round(_nu,5)==0.96432)
_gs={}
for k,d,c,Q,cls,black in _nodes:
    anc=[n for n in _nodes if n[0]<k and c%n[1]==n[2]]
    _gs[(d,c)]=sum(1.0/q for n in anc for q in n[3])+sum(1.0/q for q in Q)
_okfc=True
for N in (0.05,0.2,0.5,1.0):
    fc=[(d,c) for k,d,c,Q,cls,black in _nodes if _gs[(d,c)]>N
        and all(_gs[(n[1],n[2])]<=N for n in _nodes if n[0]<k and c%n[1]==n[2])]
    disj=all(c1%_m.gcd(d1,d2)!=c2%_m.gcd(d1,d2) for i,(d1,c1) in enumerate(fc) for (d2,c2) in fc[i+1:])
    mg=sum(1.0/d for d,c,br in _cells if sum(1.0/q for q in br)>N)
    ld=sum(1.0/cl[0] for k,d,c,Q,cls,black in _nodes if _gs[(d,c)]<=N for cl in cls)
    _okfc&=disj and abs(sum(1.0/d for d,c in fc)-mg)<1e-12 and ld<=N
ck("first crossings disjoint, measure = mu{g>N} to 1e-12, light drift <= N at four levels", _okfc)
ck("min branch survival 0.5819 > 0", round(min(_m.prod(1-1.0/q for q in br) for _,_,br in _cells),4)==0.5819)

# --- Proposition (chain constraint): exhaustive node search ---
import verify_organizable as vo
_T=[(2*q,q%2*1+2*(q%3)) for q in vo.P[2:12]]+[(10*q,4+10*(q%5)) for q in vo.P[12:20]]
_R=[(pj*q,1+pj*(q%pj)) for pj in (11,23,43) for q in vo.P[10:18] if q!=pj]
_S=[(p*q,1) for i,p in enumerate(vo.P[:6]) for q in vo.P[i+1:i+4]]
_M6=[(6*q,6) for q in vo.P[3:15]]
ck("coprime-layered: tree yes, multiples of 6 yes, protectors no, semiprimes no",
   vo.organizable(_T)[0] and vo.organizable(_M6)[0]
   and not vo.organizable(_R)[0] and not vo.organizable(_S)[0])

# --- Theorem (coprime protectors): product formula and mutation ---
import verify_coprime_protectors as vcp
_pr,_bg,_gr=vcp.instance()
_nu=1.0
for p,(c,qs) in _gr.items(): _nu*=(1-vcp.kappa(qs)/p)
ck("coprime protectors: 61 protectors, nu = 0.94479 by product", len(_pr)==61 and round(_nu,5)==0.94479)
_kap=1-_m.exp(-1.5); _pp=1.0; _prods={}
for p in vcp.primes_upto(10**6):
    if p<3: continue
    _pp*=(1-_kap/p)
    for N in (10**3,10**4,10**5,10**6):
        if p<=N: _prods[N]=_pp
ck("mutation partial products at 1e3..1e6 match", all(round(_prods[N],4)==v for N,v in [(1000, 0.2492), (10000, 0.1997), (100000, 0.168), (1000000, 0.1459)]))

# --- Theorem (overlapping protectors): exact ratios and fit ---
import verify_overlap as vov
_R=[]; _fit=[]; _nu=[]
for _N in (10**3,10**4,10**5,10**6,10**7):
    _ps=[x for x in vov.primes_upto(_N) if x>=3]
    _A,_M,_r,_n=vov.stats(_ps); _R.append(round(_r,4)); _nu.append(round(_n,4))
    _fit.append(round((_r-1)*_A/4.0,4))
ck("overlapping protectors: exact ratios [3.7216, 3.2884, 3.0348, 2.865, 2.7418]", _R==[3.7216, 3.2884, 3.0348, 2.865, 2.7418])
ck("normalised fit (R-1)A/4 falls to 1: [1.1553, 1.1345, 1.1217, 1.1129, 1.1065]",
   all(abs(a-b)<5e-4 for a,b in zip(_fit,[1.1553, 1.1345, 1.1217, 1.1129, 1.1065]))
   and all(b<a for a,b in zip(_fit,_fit[1:])) and 1.0<_fit[-1]<1.2)
ck("protector-free measure falls: [0.4811, 0.3965, 0.3392, 0.2975, 0.2658]", _nu==[0.4811, 0.3965, 0.3392, 0.2975, 0.2658])

# --- Coverage audit of the paper's catalogue ---
import io, contextlib, verify_coverage as vcov
_buf=io.StringIO()
with contextlib.redirect_stdout(_buf): _rc=vcov.main()
ck("coverage audit: all eleven catalogued systems satisfy their cited theorem's hypotheses", _rc==0)

# --- Theorem (saturation): exact reduction on finite systems ---
import verify_saturation as vsat
_b=io.StringIO()
with contextlib.redirect_stdout(_b): _rs=vsat.main()
ck("saturation reduction: contains, union, idempotent, saturated; many-protector -> 2 protectors", _rs==0)

print()
if FAILED:
    print("FAILED: %s"%", ".join(FAILED)); sys.exit(1)
print("ALL NUMERICAL CLAIMS IN THE PAPER VERIFIED.")
