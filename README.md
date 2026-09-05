# Computational artifact — *A profinite reduction of Erdős Problem 25*

This repository contains the paper and the programs behind every
numerical assertion in it.

**No theorem in the paper depends on these computations.** Each
is proved in the text; the programs exhibit the explicit
examples and counterexamples, and confirm the finite
verifications quoted in the appendix.

## Contents

- `main.tex (compiled copy: main.pdf, built with tectonic)` — the paper.
- `tools/verify_paper_claims.py` — reproduces, in one run, the values
  on which the paper's explicit assertions rest; the illustrative
  measurements quoted in remarks come from the individual probes. Several checks carry a
  control designed to fail if the property being checked were
  absent; the script exits non-zero if any check fails.
- `tools/audit_tex.py` — internal-consistency audit of the
  manuscript (reference words, dangling references, macros,
  environment and math-mode balance). It self-tests against
  seeded faults and aborts if it fails to catch them.
- `probes/` — the individual experiments, each documented in its
  own docstring, including the constructions behind the four
  refuted routes of the section of refuted routes.

## Reproducing

```
python3 tools/verify_paper_claims.py
python3 tools/audit_tex.py main.tex
```

Python 3.8 or later; standard library only, no third-party
dependencies. The claims suite runs in a few minutes and uses
well under 200 MB.

## What is verified

The exact maximum in the disjointness bound (branch-and-bound
over all 104 classes of modulus at most 14, returning exactly
1.000000); the density identities behind the incompatible-system
and coprime-divergence theorems, each against a control; the
failure of Davenport–Erdős's Lemma 1 at n = 42 for the finite
head with moduli 2..39; the cylinder reduction of the
bounded-misalignment theorem, exactly on all 83 non-empty
cylinder-class pairs; the numbers showing why the
divisibility-monotone weight does not bridge to F; and the
Heilbronn–Rohrbach enumeration over all 53,130 three-class
systems with moduli from {2,…,12}; the block measures at N = 20,
40, 80 and the transient-spike display.

Every other number in the paper — the Monte Carlo and window
measurements quoted in remarks as illustration — is produced by
the individually named program in `probes/`, run without
arguments.  The recorded output of every probe is in
`PROBE_OUTPUT/` (one file per program, from a single run of all
36 programs), so the full numerical record of the paper can be
checked against this repository without executing anything.  One
program, `verify_approx_cyl.py`, tests, on a 118-class
subsystem and moduli that are not multiples of its obstruction's
lcm, whether the unresolved share of cylinders shrinks as the
modulus grows, and records that it does not (its check line reads
FAILED). The instance lies outside the hypothesis of the paper's
remark on the relaxed decomposition (which needs M a multiple of
the obstruction's lcm and small Σ1/n_j), so the recorded failure is
uninformative rather than contradictory; the paper cites no number
from it. It is kept for the record.

## Contact

Ibrahim Mian — ibrahimnmian@gmail.com

Authors: Ibrahim Mian and Shayaan Siddique (Millennium Research). Corresponding author: Ibrahim Mian, ibrahimnmian@gmail.com.
