#!/usr/bin/env python3
"""Session 92: attempt to REFUTE Erdos 25.

The upper bound dbar_log(A) <= nu is proved, so a counterexample
needs the killed set to exceed its own measure in upper
LOGARITHMIC density. The killed set only grows, so a window's
killed density can never drop below what small-modulus classes
guarantee: liminf is pinned near 1-nu and all the room is above.

The only classes that create a TRANSIENT spike are those with
modulus comparable to X: a class with n in (X/2, X] kills
exactly one element of that window. So K such classes spike the
window by K/(X/2) while adding only K/X to the long-run drift -
the spike is about TWICE the drift.

That is the tension to measure. Survival requires the block
drifts to be summable; if the spike is bounded by twice the
drift, the spikes must vanish exactly when the survivors live,
and this refutation route is closed.

REGISTERED: REFUTATION CANDIDATE iff F oscillates with amplitude
bounded away from 0 while survivors persist. ROUTE CLOSED iff
peak/drift stays bounded and survival drives both to 0.
CONTROL: an aligned system, whose log density D-E guarantees,
must not oscillate persistently.
"""
import math


def sieve(cls, X):
    m = bytearray(X)
    for n, a in cls:
        s = a % n
        if s == 0:
            s = n
        while s < n:
            s += n
        if s < X:
            m[s::n] = b"\x01" * len(m[s::n])
    return m


def build_spike(X, scales, frac):
    """at each scale Y, insert classes of modulus in (Y/2, Y]
    aimed at DISTINCT integers of that window - the maximum
    possible transient spike"""
    cls = []
    for Y in scales:
        lo, hi = Y // 2, Y
        want = int((hi - lo) * frac)
        target = lo
        for n in range(lo + 1, hi):
            if want <= 0:
                break
            if target < n:
                target = n
            if target >= hi:
                break
            cls.append((n, target % n))
            target += 1
            want -= 1
    return sorted(set(cls))


def profile(cls, X):
    m = sieve(cls, X)
    out = []
    Y = 2048
    while Y < X:
        lo = Y // 2
        w = sum(m[lo:Y]) / float(Y - lo)
        out.append((Y, w))
        Y *= 2
    acc, F = 0.0, []
    Yi, i = 2048, 0
    for x in range(2, X):
        if m[x]:
            acc += 1.0 / x
        if i < len(out) and x == out[i][0]:
            F.append((x, acc / math.log(x)))
            i += 1
    return out, F


def main():
    X = 4 * 10 ** 6
    scales = [8192, 65536, 524288, 4194304 // 2]
    print("[spike construction] fresh classes of modulus ~Y at "
          "Y = %s, each killing one element of its own window"
          % scales)
    print()
    for frac in (1.0, 0.5, 0.2):
        cls = build_spike(X, scales, frac)
        if not cls:
            continue
        drift = sum(1.0 / n for n, _ in cls)
        w, F = profile(cls, X)
        wv = [x[1] for x in w]
        Fv = [x[1] for x in F]
        peak = max(wv)
        amp = max(Fv) - min(Fv[len(Fv) // 2:]) if len(Fv) > 3 else 0
        surv = 1 - wv[-1]
        print("  frac=%.1f  %6d classes  drift=%.4f" %
              (frac, len(cls), drift))
        print("     window killed density: %s"
              % "  ".join("%.3f" % v for v in wv))
        print("     F(X):                  %s"
              % "  ".join("%.3f" % v for v in Fv))
        print("     peak window kill %.3f, peak/drift = %.2f, "
              "final survivor density %.3f"
              % (peak, peak / drift if drift else 0, surv))
        print()

    print("[control] aligned system (D-E guarantees its log "
          "density exists) must not oscillate persistently:")
    ctrl = [(n, 0) for n in range(2, 400)]
    w, F = profile(ctrl, X)
    Fv = [x[1] for x in F]
    print("     F(X): %s" % "  ".join("%.3f" % v for v in Fv))
    tail = Fv[len(Fv) // 2:]
    print("     late-range spread %.4f -> %s"
          % (max(tail) - min(tail),
             "settles" if max(tail) - min(tail) < 0.02
             else "OSCILLATES (control broken)"))


if __name__ == "__main__":
    main()
