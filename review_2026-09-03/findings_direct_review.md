# Direct review findings (Claude Fable 5.1, 2026-09-03)

Scope: paper/main.tex, paper/supp.tex, paper/companion.tex, paper/refs.bib, read in full;
mathematics re-derived by hand and re-verified by two from-scratch scripts in this folder
(`indep_check_chsh.py`, `indep_check_companion.py`; outputs `out_indep_*.txt`).
Nothing outside this folder was modified.  Severity: **S1** = must fix before submission,
**S2** = should fix, **S3** = minor / stylistic / suggestion.

## A. Mathematical core of the main paper — verdict

Every theorem I could test holds.  Independently confirmed (see `out_indep_chsh.txt`, ALL PASS):

- Pattern census (8 even-parity patterns, four miss types, parity lemma).
- Lemma 1 (per-row TV optimum) against an LP oracle, 6900 cases, N=2..7.
- Theorem 1 + Theorem 5 curves: brute force over all 2^{4N} response vertices with a per-row
  LP (no multiset reduction, no Lemma 1) for N=2,3; exhaustive multiset formula for N=2..8 on a
  grid that includes midpoints of every cell (guards against envelope breakpoints off the grid).
- Cap theorem min_Q max_ab n_ab = ceil(N/4) for N=2..8.
- Headline: S_max(N,F*) = 2√2 exactly and S_max(N,F*−10⁻⁶) < 2√2 for N=2..8, in exact Q(√2).
- Sec. VI witness: S=2√2, E[ab]=±√2/2, every row TV = F*, pairwise TVs ∈ {F*, 2F*}.
- NS witnesses (Thm 7): the printed N≥4 closed form for N=4..13 and the printed N=3 table:
  rows valid, S=2√2, TV as stated, all four NS residuals exactly 0.
- NS N=2 flatness: LP over all 256 vertices with NS equalities, S_NS(2,F)=2 at F∈{0,.1,.25,.4,.5};
  the one-sided-NS "side finding" 2+4F also confirmed. S_NS(3,F)=2+6F confirmed by LP at 3 points.
- Thm 8 S=4 NS witnesses N=3,4,5 as printed: S=4, NS residuals 0, max row TV = ceil(N/4)/N.
- One-sided (Prop. 2) curves by brute force LP for N=2,3.
- Mutual-information numbers 0.056036 / 0.056036 / 0.054810 bits reproduced.

## B. Defects found in main.tex

**F1 (S1) Duplicate LaTeX label → wrong cross-reference in the PDF.**
`\label{eq:ns4rows}` is defined twice (line 677, the NS-witness rows for N≥4; line 838, the
S=4 class-mass table).  The reference at line 693 ("Substituting Eqs. (ns4gauge) and (ns4rows)")
resolves to the *second* definition: the compiled PDF reads "Substituting Eqs. (15) and (21)",
where (21) is the S=4 table and the intended equation is (16).  The project's "compiles clean"
criterion (zero overfull boxes, zero undefined refs) does not catch multiply-defined labels.
Fix: rename the second label (e.g. `eq:ns4caprows`).

**F2 (S2) False edge case in the one-sided remark (Sec. VII.2, "Consequences").**
"The cap to S=4 under one-sided dependence is ⌈N/2⌉/N: exactly twice the full-model cap for
N≡0,3 (mod 4), and a factor ≥3/2 (tending to 2) for N≡1,2 (mod 4)."  At N=2 the ratio is
(1/2)/(1/2)=1, not ≥3/2.  The statement is correct for N≥5.  Fix: add "N≥5" (or "N≥3").

**F3 (S2) Notation clash makes a sentence wrong as written (Sec. VII.5, "Model-class
identification").**  "Its causal structure imposes the additional factorization constraint
p(x,y|λ)=p(x|λ)p(y|λ), which our closed-form witnesses do not satisfy."  In main.tex x,y are
the *outcomes* (Sec. II), and for deterministic responses p(x,y|λ) is a point mass that trivially
factorises — so in the paper's own notation the sentence is false.  The intended statement (in
Hall–Branciard's notation, where x,y are the *settings*) is p(a,b|λ)=p(a|λ)p(b|λ); I checked
that this does fail for the witnesses (defect |p00p11−p01p10| ≈ 0.032 / 0.032 / 0.042 for the
Sec. VI, NS N=4 and NS N=3 witnesses under uniform settings).  Fix: write p(a,b|λ)=p(a|λ)p(b|λ)
and say explicitly that a,b are the settings.

**F4 (S2) Inconsistent novelty cutoff date.**  Sec. I and Sec. VII.3 say "through 2026-08-23";
the project's later literature delta check (review/…, lit_check2 Section D, memory of record)
extends the window to 2026-09-02.  Either update the manuscript date or state that the second
sweep confirmed the earlier statement.  (Also: supp.tex Sec. D says 2015–2026-08-23.)

**F5 (S3) Uniform-baseline artefacts presented as findings.**  The non-monotone cap
(1/4 at N=4, 2/5 at N=5) and the bends at F=1/N for N=5,6,7 exist only because the baseline ρ
is forced to be uniform on *all* N states (a state cannot be left unused).  With any baseline ρ
the cap is 1/4 for every N≥4 and the curve is 2+8F up to 1/4 (embed the N=4 model).  The paper
is careful to say "uniform-source" in most places, but the remark after Theorem 3 ("Note that
the cap is not 1/4 for all N≥4") and the Theorem 3 statement ("the smallest fine-tuning that
lets a local deterministic model reach S=4") read as physical statements.  A referee will call
this an artefact of the model class.  Suggest one sentence making the point explicitly.

**F6 (S3) Odd phrasing.**  Sec. VII.3: "Hall's metric is pairwise (no fixed baseline) and his
setting has no free conditional source" — Hall's model does have a settings-dependent
p(λ|a,b); the sentence is unclear.  Also "1312 candidate affine functions" (Thm 6 proof) is an
implementation detail that means nothing to a reader; either explain or drop.

**F7 (S3) Trivial-lemma attribution.**  The "PR-box uniqueness" lemma is standard (the PR box
is the unique NS box saturating S=4); it should be attributed (Popescu–Rohrlich 1994 is not
even in refs.bib) rather than presented as new.

**F8 (S3) Placeholders.**  Author block, acknowledgments, companion citation, code URL
(supp.tex) are all TODO — expected, but the companion paper is cited as a real reference for
several results (Theorem 8 companion link, K/m_min conjecture), so the two papers must be
submitted together or the main paper should not lean on it.

## C. Companion paper (companion.tex) — findings

**F9 (S1) The CGLMP "quantum value" used is not the quantum maximum.**  The companion takes
T=(12+8√3)/9≈2.87293 (CGLMP's maximally-entangled-state value) as the target and presents
F*=√3/9−1/12≈0.109117 as the headline.  A larger qutrit value is well known: Acín, Durt, Gisin,
Latorre, PRA 65, 052325 (2002) obtain I₃ = 1+√(11/3) ≈ 2.9149 with the non-maximally entangled
state (|00⟩+γ|11⟩+|22⟩)/√(2+γ²), γ=(√11−√3)/2, and the same CGLMP measurements; this is also the
NPA/Tsirelson value for I₃ found numerically in later work.  My script reproduces both numbers
(`out_indep_companion.txt`).  Consequences: (i) the abstract's "the achieved quantum value" is
misleading — it is the *maximally-entangled* value, not the achieved quantum maximum; (ii) with
the correct target, F*(N≥4) = (√(11/3)−1)/8 ≈ 0.11436 (still on the slope-8 segment, so the
formula F*=(T−2)/(2 min(N,4)) survives, but the quoted number changes); (iii) the "prior guess
≈3.34 not confirmed" remark in Sec. VIII should be replaced by the correct citation.

**F10 (S2) False general statement in Sec. II.C.**  "listing V(k/N) for k=0,…,N gives the
complete curve."  The envelope over multisets can bend *between* grid points: for Mermin n=5,
N=2 the exact values are V(0)=4, V(1/4)=8, V(1/2)=16 (my script), so the curve is 4+16F on
[0,1/4] and 32F on [1/4,1/2], not the chord.  The paper's own Sec. V text acknowledges this
bend, contradicting Sec. II.C.  The PM-square convexity argument (Sec. VI) is the correct
version and is valid; Sec. II.C should be weakened to "when the round-robin curve matches the
envelope at all grid points".

**F11 (S3) Theorem (Mermin n=5) is weaker than what is provable.**  Equality
F*=⌈3N/8⌉/N is claimed for N≡0,1,2 (mod 16).  Since the exact values at N=3,4,6 also equal
⌈3N/8⌉/N, "m copies of all 16 min-miss patterns + the optimal size-r multiset" gives equality
for r∈{0,1,2,3,4,6} as well; only r=5 and r∈{7..15} are open.  Not an error, but the sandwich
upper bound (6m+r)/(16m+r) is needlessly loose.

**F12 (S3) Mermin n=3 is literally the CHSH problem.**  The pattern structure (8 even-parity
patterns, values ±2, four miss types) is isomorphic to CHSH, so Theorem (Mermin n=3) is a
re-labelling of the main paper's cap theorem.  Fine to include, but say so.

## D. Verified numbers in the companion (my script)

Mermin n=3: 8 patterns, B=2, m_min=1, 4 near-miss.  Mermin n=5: 32 patterns, B=4, m_min=6,
16 min-miss with incidence 6, A_maxB=(8,12,16,16,16) for N=2..6, F*=(1/2,2/3,1/2,3/5,1/2),
V(k/N) tables as printed.  GHZ₃/GHZ₅ give 4/16 with perfect correlators.  CGLMP: 27 patterns,
local bound 2, four near-miss patterns as listed, generalized row lemma vs LP, F*(T;N) for
N=2,3,4.  PM square: line products (+,+,+,−,+,+), 32 odd miss-sets, caps ⌈N/6⌉/N and the exact
curve table for N=2..10 (DP over achievable w-vectors, checked on grid + midpoints).
(See `out_indep_companion.txt` for the final status of each item.)
