#!/usr/bin/env python3
"""Systematic consistency audit for the Erdos 25 paper.

Written because three edits failed SILENTLY across recent
sessions, one leaving the vanishing-tails proof stating
Davenport-Erdos for D^{(N)} while using D_N two lines later -
an inconsistency that survived two adversarial review rounds.
Verifying each new edit catches future failures; this finds
damage already done.

Checks:
 1. REFERENCE WORDS. Record the environment each \\label was
    declared in, then check every "Word~\\ref{label}" names that
    environment. Downgrading a theorem to a proposition and back
    is exactly how stale words are left behind.
 2. UNDEFINED MACROS. A half-applied edit can leave a command
    whose \\newcommand was removed.
 3. USE BEFORE DEFINITION for the paper's core symbols.
 4. Environment balance and dangling references.

The checker is SELF-TESTED against seeded faults: a checker that
reports "clean" because it is broken is worse than none.
"""
import re
import sys

ENVW = {"theorem": "Theorem", "lemma": "Lemma",
        "proposition": "Proposition", "remark": "Remark",
        "definition": "Definition", "corollary": "Corollary",
        "equation": "equation"}


def label_envs(t):
    """map label -> environment it was declared in"""
    out = {}
    for m in re.finditer(r"\\begin\{(\w+)\}(\[[^\]]*\])?"
                         r"\\label\{([^}]+)\}", t):
        out[m.group(3)] = m.group(1)
    # equations labelled inside display math
    for m in re.finditer(r"\\begin\{equation\}\\label\{([^}]+)\}", t):
        out[m.group(1)] = "equation"
    return out


def check_ref_words(t, envs):
    bad = []
    # connectives follow a plural type word ("Lemmas~X and~Y")
    SKIP = {"and", "or", "of", "in", "by", "from", "see", "to",
            "with", "via"}
    for m in re.finditer(r"(\w+)~\\ref\{([^}]+)\}", t):
        word, lab = m.group(1), m.group(2)
        if lab not in envs or word in SKIP:
            continue
        want = ENVW.get(envs[lab])
        # irregular plural: "Corollaries", not "Corollarys"
        PLUR = {"Corollary": "Corollaries"}
        plurals = {want + "s", PLUR.get(want, want + "s")}
        if want and word not in ({want} | plurals):
            bad.append((word, lab, want))
    return bad


def check_macros(t):
    defined = set(re.findall(r"\\newcommand\{\\(\w+)\}", t))
    known = set("""begin end label ref eqref cite section subsection
      documentclass usepackage newtheorem theoremstyle title
      author email maketitle emph textup textbf item bibitem
      frac sum prod lim limsup liminf sup inf max min log exp
      le ge ne in subseteq subset cup cap setminus sqcup mid nmid
      equiv pmod bmod cdots ldots dots to infty alpha beta delta
      varepsilon varphi zeta nu mu sigma Lambda quad qquad
      text mathbb mathrm widehat varprojlim colon left right
      bigl bigr Bigl Bigr bigcup bigcap displaystyle
      thebibliography medskip noindent varnothing
      H o u i c s appendix nonumber notag qed
      tfrac binom over choose hspace vspace par
      forall exists neq leq geq approx sim ll gg
      Longrightarrow longrightarrow cdot cong delta iff notin
      gcd lcm gtrsim lesssim operatorname mathcal
      gamma int lceil lfloor rceil rfloor
      overline underline downarrow emptyset times prod
      date normalsize texttt url Sigma asymp ell eta kappa mapsto
      ni not pi rho sqrt supseteq
      """.split())
    # letters only: \w would slice \sum_{i>N} into "sum_"
    used = set(re.findall(r"\\([A-Za-z]+)", t))
    return sorted(used - defined - known)


def check_symbol_order(t, symbols):
    body = t[t.index("\\section{Introduction}"):]
    bad = []
    for sym, definer in symbols:
        first = body.find(sym)
        dpos = body.find(definer)
        if first >= 0 and dpos >= 0 and first < dpos:
            bad.append((sym, first, dpos))
    return bad


def check_math_mode(t):
    """balance checks a compiler would catch; we have no LaTeX
    toolchain here, so do them by hand"""
    msgs = []
    body = t.replace("\\$", "")
    if body.count("$") % 2:
        msgs.append("odd number of $ (unbalanced inline math)")
    # a line break with an optional length, \\[2pt], is not display math
    opens = len(re.findall(r"(?<!\\)\\\[", body))
    closes = len(re.findall(r"(?<!\\)\\\]", body))
    if opens != closes:
        msgs.append("unbalanced display math \\[ \\]: %d vs %d" % (opens, closes))
    if body.count("\\bigl") != body.count("\\bigr"):
        msgs.append("unbalanced \\bigl/\\bigr: %d vs %d"
                    % (body.count("\\bigl"), body.count("\\bigr")))
    if body.count("\\Bigl") != body.count("\\Bigr"):
        msgs.append("unbalanced \\Bigl/\\Bigr")
    if body.count("\\left") != body.count("\\right"):
        msgs.append("unbalanced \\left/\\right")
    return msgs


def check_structure(t):
    msgs = []
    for env in ("document", "abstract", "theorem", "proof",
                "definition", "lemma", "proposition", "remark",
                "thebibliography"):
        o = len(re.findall(r"\\begin\{" + env + r"\}", t))
        c = len(re.findall(r"\\end\{" + env + r"\}", t))
        if o != c:
            msgs.append("unbalanced %s: %d/%d" % (env, o, c))
    labs = set(re.findall(r"\\label\{([^}]+)\}", t))
    refs = set(re.findall(r"\\ref\{([^}]+)\}", t))
    refs |= set(re.findall(r"\\eqref\{([^}]+)\}", t))
    if refs - labs:
        msgs.append("dangling refs: %s" % sorted(refs - labs))
    msgs += check_math_mode(t)
    return msgs


def run(t, quiet=False):
    envs = label_envs(t)
    rw = check_ref_words(t, envs)
    mc = check_macros(t)
    so = check_symbol_order(t, [
        ("D_c", "D_c="), ("u_i", "n=uw"), ("K_i", "K_i="),
    ])
    st = check_structure(t)
    if not quiet:
        print("  labels found: %d" % len(envs))
        print("  stale reference words: %s"
              % ("; ".join("'%s~\\ref{%s}' should be '%s'"
                           % b for b in rw) if rw else "none"))
        print("  unknown macros: %s" % (mc if mc else "none"))
        print("  symbols used before definition: %s"
              % (so if so else "none"))
        print("  structural: %s" % ("; ".join(st) if st else "ok"))
    return rw, mc, so, st


def self_test(t):
    """seed faults; the checker MUST catch each one"""
    print("[self-test] seeding faults the checker must catch")
    ok = True
    # fault 1: wrong reference word
    lab = re.search(r"\\begin\{theorem\}\\label\{([^}]+)\}", t)
    if lab:
        k = lab.group(1)
        bad = t.replace("Theorem~\\ref{%s}" % k,
                        "Proposition~\\ref{%s}" % k, 1)
        rw, _, _, _ = run(bad, quiet=True)
        hit = any(b[1] == k for b in rw)
        print("   wrong reference word detected: %s" % hit)
        ok &= hit
    # fault 2: dangling reference
    bad = t.replace("\\label{lem:tail}", "\\label{lem:tail_X}", 1)
    _, _, _, st = run(bad, quiet=True)
    hit = any("dangling" in m for m in st)
    print("   dangling reference detected:    %s" % hit)
    ok &= hit
    # fault 3: unbalanced environment
    bad = t.replace("\\end{proof}", "", 1)
    _, _, _, st = run(bad, quiet=True)
    hit = any("unbalanced" in m for m in st)
    print("   unbalanced environment detected: %s" % hit)
    ok &= hit
    return ok


def main():
    t = open(sys.argv[1]).read()
    if not self_test(t):
        print("\n[ABORT] the checker failed its own self-test; "
              "its clean verdict would mean nothing.")
        sys.exit(2)
    print("\n[audit]")
    rw, mc, so, st = run(t)
    problems = bool(rw) or bool(so) or bool(st)
    print()
    if problems:
        print("[RESULT] issues found above - each must be fixed "
              "or explained.")
        sys.exit(1)
    print("[RESULT] no inconsistencies of the kinds a silently "
          "failed edit produces.")


if __name__ == "__main__":
    main()
