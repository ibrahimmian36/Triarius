#!/usr/bin/env python3
"""Session 101: the positive-harmonic-density regime.

c = limsup D(Y)/log Y.  The spike bound of session 97 is vacuous
when c > 0, so those systems need separate treatment.

TWO FACTS FOUND BEFORE CODING.
 * The doubled construction, the campaign's hardest example, has
   c = 0 (its drift grows like log log, not log). So the hard
   regime is NOT where that example lives.
 * c > 0 blocks the tools: the evens contain NO pairwise-coprime
   pair at all, so Theorem D is unavailable there by
   construction.

THE DICHOTOMY TO TEST. The evens at residue 0 have c = 1/2 and
survive with density 1/2 - but they are ALIGNED, already closed
by translation. Can a dense modulus set with NON-ALIGNED
residues also survive, or does density plus misalignment force
death?  Heuristically survival needs x to dodge one residue
class per modulus up to x, of probability prod(1-1/n) ~ 1/x,
suggesting death - but that heuristic has been wrong before, so
it is measured.

Class counts are CAPPED: last session overran to 668 MB because
a control had 2.75M classes, not because X was large.
"""
import math
import random
import sys

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)

# CAP FIX: capping the class count at 60000 stopped the moduli
# at 120000 while the window ran to 2e6, so the system was
# FINITE and its survivors plateaued (0.0052, 0.0166) - an
# artifact of the cap, not a real limit. The moduli must reach
# through the window. Classes are now STREAMED rather than
# stored, so memory stays bounded without capping the count.


def survivors(cls, X):
    m = bytearray(X)
    for n, a in cls:
        s = a % n
        if s == 0: s = n
        while s < n: s += n
        if s < X: m[s::n] = b"\x01"*len(m[s::n])
    out = []
    Y = 8192
    while Y <= X:
        out.append((Y, 1 - sum(m[Y//2:Y])/float(Y-Y//2)))
        Y *= 2
    return out


def harmonic_c(mods, Y):
    return sum(1.0/n for n in mods if n <= Y)/math.log(Y)


def main():
    X = 2*10**6
    rng = random.Random(11)
    fams = {
      "all integers":   range(2, X//2),
      "evens":          range(2, X//2, 2),
      "multiples of 3": range(3, X//2, 3),
    }
    print("dense modulus families (c = harmonic density):")
    for name, mods in fams.items():
        print("   %-16s %8d moduli, c ~ %.3f"
              % (name, len(mods), harmonic_c(mods, X)))
    print()
    print("1. ALIGNED controls - these must SURVIVE at their "
          "known densities")
    ctrl_ok = True
    for name, mods, pred in (("evens (res 0)",
                              fams["evens"], 0.5),
                             ("multiples of 3 (res 0)",
                              fams["multiples of 3"], 2/3.0)):
        s = survivors(((n, 0) for n in mods), X)
        got = s[-1][1]
        print("   %-24s survivor density %.4f (predicted %.4f)"
              % (name, got, pred))
        if abs(got-pred) > 0.02: ctrl_ok = False
    ck("aligned controls survive at the arithmetic prediction",
       ctrl_ok)

    print()
    print("2. the same families with NON-ALIGNED residues")
    verdicts = {}
    for name, mods in fams.items():
        s = survivors(((n, rng.randrange(n)) for n in mods), X)
        seq = [v for _, v in s]
        falling = all(b <= a + 1e-9 for a, b in zip(seq, seq[1:]))
        verdicts[name] = (seq, falling)
        print("   %-16s survivor density: %s   %s"
              % (name, "  ".join("%.4f" % v for v in seq[-6:]),
                 "falling" if falling else "NOT monotone"))
    print()
    dying = [n for n, (seq, f) in verdicts.items()
             if seq[-1] < 0.02]
    living = [n for n, (seq, f) in verdicts.items()
              if seq[-1] >= 0.02]
    if living:
        print("[verdict] SURVIVAL in %s (final densities %s). A "
              "dense NON-aligned system holding positive density "
              "is a genuinely new hard case and matters more "
              "than the negative would have."
              % (", ".join(living),
                 ", ".join("%.4f" % verdicts[n][0][-1]
                           for n in living)))
    else:
        print("[verdict] DEATH in every dense non-aligned family "
              "tested (%s). Density plus misalignment appears to "
              "force the survivors to 0, which would close the "
              "c > 0 regime by triviality. MEASURED, not proved."
              % ", ".join(dying))
    if FAIL:
        sys.exit(1)


if __name__ == "__main__":
    main()
