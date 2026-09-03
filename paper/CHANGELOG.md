# What the review changed in the Qwen-written papers

The manuscripts in `original_qwen_2026-09-02/` are exactly what the local Qwen3.8-27B model
produced.  The corrected manuscripts in this directory differ from them in the ways listed below.
Everything here was found and applied by Claude Fable 5.1 (with Claude Opus 5 sub-agents for the
literature and certificate audits) on 2026-09-03.  Only changes that bear on the question this
repository is about — what a small local model can and cannot do in a research project of this
kind — are listed; author fields, layout and other housekeeping are omitted.

## 1. Mathematics: confirmed, not changed

No theorem, proof, number, curve, witness or certificate was found wrong.  The review re-derived
the results by hand and by independently written exact-arithmetic and linear-programming checks
(`../review_2026-09-03/`), re-ran all 26 project certificates, and confirmed every quoted number
against its log.  Three statements *about* the mathematics did need correction:

- **Main paper, Sec. VII.2.** "The one-sided cap is a factor ≥ 3/2 above the full cap for
  N ≡ 1, 2 (mod 4)" fails at N = 2, where the ratio is 1.  Restricted to N ≥ 5.
- **Companion, Sec. II.** "Listing V(k/N) for k = 0,…,N gives the complete curve" is false in
  general: the upper envelope over multisets can bend between grid points, and the paper's own
  Mermin n = 5, N = 2 curve does so (V(1/4) = 8, not the chord value 10).  Replaced by the correct
  convexity statement.
- **Companion, Mermin n = 5 theorem.** Equality F* = ⌈3N/8⌉/N was claimed for residues 0, 1, 2
  (mod 16) only; the paper's own exact values at N = 3, 4, 6 make it provable for residues
  0, 1, 2, 3, 4, 6.  Strengthened, with the open residues corrected to {5, 7, …, 15}.

## 2. A wrong physical target

- **Companion, CGLMP d = 3.** The "quantum value" used as the target, T = (12+8√3)/9 ≈ 2.8729, is
  the value for the *maximally entangled* two-qutrit state, not the quantum maximum.  Acín, Durt,
  Gisin and Latorre (PRA 65, 052325, 2002) obtained 1+√(11/3) ≈ 2.9149 with a non-maximally
  entangled state, and the Navascués–Pironio–Acín hierarchy certifies that value as the maximum to
  numerical precision.  The model did not know this: an earlier draft even recorded an unconfirmed
  guess "≈ 3.34" for the maximum.  Both targets are now treated; the closed form is unchanged, and
  the N ≥ 4 headline number at the true maximum is 0.114357 instead of 0.109117.

## 3. Bibliography: four bad entries, two of them marked "verified"

- **Tsirelson.**  Cited as Lett. Math. Phys. 6, 293–306 (1974).  The paper is Cirel'son,
  Lett. Math. Phys. 4, 93–100 (1980).  Volume, pages and year were all wrong, in the reference
  used for the Tsirelson bound in the opening paragraph.
- **Shimony.**  Cited as "Bell's theorem and quantum field theory", eds. Mittelstaedt & Michel,
  Springer, pp. 327–340, with an in-file comment claiming an "exact match" against a source that
  in fact gives a different title, editors, publisher and pages.  Rebuilt as "Controllable and
  uncontrollable non-locality", Proc. ISQM Tokyo 1983, ed. Kamefuchi et al., Physical Society of
  Japan (1984), pp. 225–230.
- **Vieira, Ramanathan & Cabello.**  DOI had transposed digits and did not resolve.
- **Jarrett.**  Wrong end page (569–579 → 569–589).
- **Companion.**  The Peres–Mermin square was cited to Mermin's GHZ paper (PRL 65, 1838) instead
  of Mermin PRL 65, 3373 and Peres 1990.

Every other entry, and every numerical claim attributed to a cited work (Hall, Hall–Branciard,
Takakura et al., Pal et al., Alai, CGLMP, Garza–Hance), checked out against the primary source.

## 4. Literature awareness

- **Missing central references.**  Neither paper cited Barrett & Gisin (PRL 106, 100406, 2011),
  Pütz et al. (PRL 113, 190402, 2014), Pütz & Gisin (NJP 18, 055006, 2016) or Koh et al.
  (PRL 109, 160404, 2012) — the prior results on exactly the question asked, and the ones showing
  that the answer depends strongly on the metric.  Added to the introduction with a new paragraph,
  "Metric dependence of the answer", framing the result as one exactly solved point in a
  metric-dependent landscape rather than a metric-independent price of the loophole.
- **Popescu–Rohrlich** was not cited for the PR box, whose uniqueness was presented as a lemma;
  now attributed as a standard fact.
- **The headline itself has a near-equivalent in the missed literature.**  A follow-up check
  (`../review_2026-09-03/puetz_correspondence.py`) found that Pütz et al. 2014, Eq. (11) / Pütz &
  Gisin 2016, Eq. (31) — CHSH ≤ 4(1−2ℓ) for measurement-dependent local models with no-signalling
  and uniform inputs, where ℓ is the lower bound on P(x,y|λ), with quantum violation only for
  ℓ > (2−√2)/4 — becomes exactly S ≤ 2+8F with F* = (√2−1)/4 under ℓ = 1/4 − F, and that their
  extremal construction *is* the paper's round-robin N=4 model.  The metrics are not equivalent in
  general, so the TV bound for arbitrary sources and N, the finite-N curves and caps, the N=2
  no-signalling impossibility, the S=4 cap and the extensions remain new; but the four-state number
  is a re-parametrisation of a 2014 result.  The paper now states this in the abstract,
  introduction and Sec. VII.3, and its novelty claims were weakened accordingly.  For the
  experiment this is the most telling defect: the model reproduced a published result in a new
  variable without recognising or citing it.

## 5. Statements false or overstated as written

- **Main paper, Sec. VII.5.**  Hall–Branciard's causal-structure constraint was transcribed in
  their notation, p(x,y|λ) = p(x|λ)p(y|λ), in which x, y are settings; in this paper x, y are
  outcomes, for which the factorisation holds trivially, so the clause "which our witnesses do
  not satisfy" was false as written.  Rewritten in terms of the settings.
- **Alai's measure** was labelled "Hall-type variational M" in Table I but filed as Hall's M/2 in
  the text — a factor-2 inconsistency; and his work was described as "LP with Q(√2) certificates"
  when the exact certificate covers a single point.  Both corrected.
- **Vieira et al.** were grouped under "experimental literature"; it is a theoretical exclusion
  result, and one directly adjacent to this paper's topic.  Now described and engaged accurately.
- **"No floating point"** (supplement) was a blanket claim; the LP cross-checks, the
  mutual-information numerics and the Garza–Hance evaluation are floating point.  Qualified.
- **Companion, CGLMP census** mentioned only the four near-miss patterns; six further patterns of
  weight type (+1,+1,0,0) also attain the local bound.  Stated.

## 6. Bookkeeping that affected the record

- A duplicated equation label made the proof of Theorem 7 cite the wrong equation in the PDF.
- The literature-sweep cutoff date in the manuscripts (2026-08-23) had not been advanced after
  the project's own later delta check (2026-09-02).
- An existing exhaustive certificate for N = 9,…,14 (`curves_highN.py`) was not cited; both
  papers said N ≥ 9 was "left to numerics".  It is now cited, and was independently re-derived.
- The abstracts listed every result at equal weight; they were rewritten to state the
  contribution and its scope first.

## 7. Added, not corrected

- A "Note on provenance and corrections" at the end of each paper.
- The former Supplemental Material folded into the main paper as Appendices A–D (for arXiv).
- The positioning paragraph and references of item 4; the CGLMP maximum and its two references
  of item 2; the uniform-baseline remark explaining why the finite-N caps are non-monotone.
