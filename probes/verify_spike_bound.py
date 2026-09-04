#!/usr/bin/env python3
"""Session 97: a quantitative bound on the spike part of the
excursion set.

FACTS. (i) A class of modulus n in (X/2, X] kills AT MOST ONE
element of the window [X/2, X], since x >= n and x <= X < 2n.
(ii) The moduli are distinct, so with N counting those classes,
N/X <= sum_{n in (X/2,X]} 1/n: the spike contribution 2N/X is at
most twice the window's own drift.

THEOREM. Summing window drifts up to Y gives D(Y), and there are
log Y / log 2 windows, so the mean spike contribution is
2 log2 * D(Y)/log Y -> 2 log2 * c, where c = limsup D(Y)/log Y
is the harmonic density of the modulus set. By Markov the
logarithmic density of scales whose spike exceeds delta is at
most (2 log 2) c / delta - zero when c = 0.

SCOPE, stated honestly: this bounds the SPIKE term only, the
contribution of moduli comparable to X. Classes of intermediate
modulus are NOT controlled this way; a union bound over them
gives about c log X, and refining by the old survivor density
gives nu_T c log X, both useless.

CHECKS. (1) the one-element fact, by direct enumeration;
(2) spike <= 2 * window drift, every window; (3) a control with
c = 0 must show vanishing spikes, and one with large c must show
them - else the bound is untested.
"""
import math
import sys

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def check_one_element(X=200000):
    """a class of modulus in (X/2,X] kills at most one element
    of [X/2,X] - enumerated, not assumed"""
    worst = 0
    tested = 0
    for n in range(X//2+1, X, 997):
        for a in (0, 1, n//3, n-1):
            s = a % n
            if s == 0: s = n
            while s < n: s += n
            cnt = len([x for x in range(s, X, n) if x >= X//2])
            worst = max(worst, cnt)
            tested += 1
    ck("class of modulus in (X/2,X] kills at most one window "
       "element", worst <= 1 and tested > 50,
       "max %d over %d classes" % (worst, tested))
    return worst <= 1


def windows(cls, X):
    """per-window: spike contribution and window drift"""
    out = []
    Y = 4096
    while Y <= X:
        lo, hi = Y//2, Y
        big = [(n, a) for n, a in cls if lo < n <= hi]
        killed = set()
        for n, a in big:
            s = a % n
            if s == 0: s = n
            while s < n: s += n
            for x in range(s, hi, n):
                if x >= lo: killed.add(x)
        spike = len(killed)/float(hi-lo)
        drift = sum(1.0/n for n, _ in big)
        out.append((Y, spike, drift))
        Y *= 2
    return out


def main():
    X = 2**21
    print("1. the one-element fact")
    check_one_element()
    print()
    print("2. spike <= 2 * window drift, and the harmonic "
          "density c")
    systems = {
      "all moduli (c ~ 1)": [(n, n//2) for n in range(2, X)],
      "even moduli (c ~ 1/2)": [(n, 0) for n in range(2, X, 2)],
      "squares (c = 0)": [(k*k, 1 % (k*k))
                          for k in range(2, int(X**0.5))],
    }
    for name, cls in systems.items():
        w = windows(cls, X)
        viol = [(Y, s, d) for Y, s, d in w if s > 2*d + 1e-9]
        D = sum(1.0/n for n, _ in cls if n <= X)
        c = D/math.log(X)
        print("  %-22s c ~ %.3f   spikes: %s"
              % (name, c,
                 "  ".join("%.3f" % s for _, s, _ in w[-5:])))
        ck("  spike <= 2*(window drift) in every window [%s]"
           % name, not viol,
           "%d violations" % len(viol))
    print()
    print("3. the bound must discriminate: c=0 gives vanishing "
          "spikes, large c gives large ones")
    sp0 = max(s for _, s, _ in windows(systems["squares (c = 0)"], X))
    sp1 = max(s for _, s, _ in windows(systems["all moduli (c ~ 1)"], X))
    ck("c=0 system has negligible spikes", sp0 < 0.02,
       "max spike %.4f" % sp0)
    ck("c~1 system has large spikes (so the test is not "
       "vacuous)", sp1 > 0.3, "max spike %.4f" % sp1)
    print()
    if FAIL:
        print("FAILED: %s" % ", ".join(FAIL)); sys.exit(1)
    print("[verdict] spike term bounded by twice the window "
          "drift; Markov then bounds the log-density of "
          "spike-excursions by (2 log 2) c / delta.")


if __name__ == "__main__":
    main()
