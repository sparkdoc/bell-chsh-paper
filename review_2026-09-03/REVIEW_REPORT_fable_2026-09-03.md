> **Note for readers of this repository.**  This is the review report as written in the working
> repository.  Paths it mentions map as follows: `paper_fixed/` is `paper/` here; `paper/` (the
> unreviewed originals) is `paper/original_qwen_2026-09-02/`; `patch_paper*.py` are the scripted
> edits that produced the corrected sources (not shipped; the resulting changes are summarised in
> `paper/CHANGELOG.md`); `rerun/` and `build/` evidence directories are not shipped.

# Independent critical review — Bell/CHSH fine-tuning project (2026-09-03)

Reviewer: Claude Fable 5.1 (main session) with two Claude Opus 5 sub-agents
(literature/bibliography; certificates/records/compile).  All prior project content was
produced by Qwen3.8-27B; this review shares no code with the project's own certificates.
Nothing outside `review/fable_2026-09-03/` was modified (`git status` unchanged except the
pre-existing deletion of `HANDOFF_tier3_fold.md`).

## Files in this folder

| file | content |
|---|---|
| `REVIEW_REPORT_fable_2026-09-03.md` | this consolidated report |
| `findings_direct_review.md` | my line-by-line findings on main.tex / supp.tex / companion.tex |
| `indep_check_chsh.py` → `out_indep_chsh.txt` | from-scratch verification of the main paper (ALL PASS) |
| `indep_check_companion.py` → `out_indep_companion.txt` | from-scratch verification of the companion (ALL PASS after fixing two convention slips in *my* script; first run kept as `out_indep_companion_run1.txt`) |
| `cglmp_quantum.py` → `out_cglmp_quantum.txt` | CGLMP quantum value: maximally entangled vs ADGL state |
| `lit_A_bibliography.md`, `lit_B_claims.md` | Opus agent: every refs.bib entry; ~45 claims about cited works |
| `cert_audit.md`, `docs_consistency.md`, `compile_check.md` | Opus agent: 15 certificate scripts re-run and audited; PROOFS/RESULTS/REVIEW_REPORT/README vs papers; compile diagnostics |
| `rerun/`, `build/` | evidence: re-run logs, regenerated figure, scratch PDF builds |

## 1. Bottom line

**The mathematics of the main paper is sound.**  Every theorem, lemma, curve, cap, witness and
NS statement I could test was re-derived by hand and reproduced by independent computation
(LP oracles for Lemma 1 and the one-sided lemma, brute force over all response vertices for
N=2,3, exhaustive multiset enumeration for N=2..8 including off-grid midpoints, exact Q(√2)
checks of every printed witness, LP-over-all-vertices for the NS N=2 flatness).  The project's
own 15 certificate scripts all re-ran to exit 0 and reproduce their logs; all 38 numbers the
manuscripts quote from those logs match.  All headline numbers attributed to Hall,
Hall–Branciard, Takakura, Pal, Alai, CGLMP and Garza–Hance were confirmed against primary
sources, several by recomputation.

**The papers are nonetheless not submission-ready.**  There is one substantive scientific
problem in the companion paper (wrong CGLMP quantum target), one wrong cross-reference inside
a proof in the main paper, four bibliography entries with real errors (two of them carrying
"verified" comments that are demonstrably false), a false sentence caused by a notation
clash, and a number of smaller record-keeping inconsistencies.  Hallucination pattern
observed: the model's *math* is reliable; its *bibliographic provenance claims* are not.

## 2. Must fix before submission (S1)

| id | where | problem | fix |
|---|---|---|---|
| **F9** | companion.tex abstract, Sec. III, Sec. VIII | The CGLMP "quantum value" T=(12+8√3)/9≈2.87293 is the *maximally-entangled-state* value, not the quantum maximum.  Acín–Durt–Gisin–Latorre, PRA 65, 052325 (2002) obtain I₃=1+√(11/3)≈2.914854 with the non-maximally entangled state (|00⟩+γ|11⟩+|22⟩)/√(2+γ²), γ=(√11−√3)/2 — reproduced numerically here (`out_cglmp_quantum.txt`; my scan over γ finds no larger value in that family).  The abstract's phrase "the achieved quantum value" is misleading and the Sec. VIII remark about an unconfirmed "≈3.34" shows the authors did not know the ADGL result. | Cite ADGL; either adopt T'=1+√(11/3) as the target (F*(N≥4)=(√(11/3)−1)/8≈0.114357; the formula F*=(T−2)/(2min(N,4)) still applies since T' lies on the slope-8 segment) or state explicitly that T is the maximally-entangled value and give both numbers. |
| **F1 / C-1** | main.tex:677 and :838 | `\label{eq:ns4rows}` defined twice.  The reference at :693 ("Substituting Eqs. (15) and (21)") resolves to the S=4 table (21) instead of the NS-witness rows (16) — a wrong cross-reference inside the proof of Theorem 7. | Rename the second label. |
| **A-3** | refs.bib `tsirelson1974` | Volume, pages, year and key all wrong.  Correct: B. S. Cirel'son, Lett. Math. Phys. **4**, 93–100 (**1980**), doi 10.1007/BF00417500. | Replace entry; cited in the paper's opening paragraph. |
| **A-5** | refs.bib `shimony1984` | Entry is garbled: title "Bell's theorem and quantum field theory", eds. Mittelstaedt & Michel, Springer, pp. 327–340 — none of these match the source the in-file comment says it was verified against (Todorov quant-ph/9701024 ref. [11], read from the local PDF).  Correct: A. Shimony, "Controllable and uncontrollable non-locality", in *Foundations of Quantum Mechanics in the Light of New Technology*, ed. S. Kamefuchi et al. (Physical Society of Japan, Tokyo, 1984), pp. 225–230. | Rebuild entry; correct the false "C1 RESOLVED" comment and the memory line. |
| **A-12** | refs.bib `vieira2025` | DOI does not resolve (digits transposed).  Correct: 10.1038/s41467-025-59247-7. | Fix DOI. |

## 3. Should fix (S2)

Main paper text
- **F3 / B2-c** (main.tex Sec. VII.5): "p(x,y|λ)=p(x|λ)p(y|λ), which our witnesses do not satisfy" is false in the paper's own notation (x,y are outcomes; deterministic models trivially factorise them).  Intended: p(a,b|λ)=p(a|λ)p(b|λ) with a,b the settings — I verified the witnesses do violate that (defect ≈0.03–0.04).
- **F2** (Sec. VII.2): one-sided/full cap ratio "≥3/2 for N≡1,2 (mod 4)" fails at N=2 (ratio 1).  Add N≥5.
- **B6-e / B6-d** (Table I, Sec. VII.3): Alai's M is Hall's M/2 (fixed by his η=1 edge), but Table I labels it "Hall-type variational M" while the disambiguation files it as M/2 — factor-2 inconsistency; and "LP with Q(√2) certificates" overstates Alai (exact certificate at one point η=9/10 only).
- **B14-d**: Vieira–Ramanathan–Cabello is a theory paper, not "experimental literature"; its result (partial MI relaxations are excludable with many settings/high dimension) is adjacent to this paper's topic and deserves a sentence.
- **F4 / B17 / D-3**: novelty cutoff.  Manuscripts, supp and RESULTS say 2026-08-23; the 2026-09-02 delta check exists only in its log; the memory says 09-02.  Pick one and apply everywhere.
- **D-10**: main.tex:609,733 cite "Supplemental Material, Sec. C" but the supplement numbers its sections I–IV (the verification record is Sec. III).
- **F6**: "his setting has no free conditional source" (Hall) is unclear/incorrect; "1312 candidate affine functions" is an unexplained implementation detail.
- **A-4**: Jarrett pages 569–**589**.  **C-4**: `[Author name]` renders as "A. name]".

Companion paper text
- **F10**: Sec. II.C "listing V(k/N) … gives the complete curve" is false in general — for Mermin n=5, N=2 the envelope bends at F=1/4 (V(1/4)=8, not the chord value 10); the paper's own Sec. V describes that bend.  The PM-square convexity argument (Sec. VI) is the correct form.
- **F13**: CGLMP census.  Ten of the 27 patterns attain the local bound 2 (four near-miss of type (+,+,+,−) and six of type (+,+,0,0)); the text mentions only the four.  Results unaffected (Prop. cap correctly counts w<+1), but the "all-but-one structure" wording should acknowledge the 0-weight patterns.
- **B16-g**: `mermin1990` (PRL 65, 1838, the GHZ paper) is cited for the Peres–Mermin square; canonical is Mermin PRL 65, 3373 (1990) (and Peres 1990).
- **B7-c**: "each odd-n value repeated at the following even size" — Alai says "through n=11".
- **C-7**: `\lfloor N/6\rceil` typo (companion.tex:777); "F ≥ max_c w_c(Q)" missing "/N" (:843).
- **C-2**: companion uses `\bibliographystyle{apsrev}` vs main's `apsrev4-2`.

Certificates / figure / records (from `cert_audit.md`, `docs_consistency.md`)
- **A-17**: in `paper/fig_curves.pdf` the N=4 curve is invisible (N=8 has the identical curve and overplots it) while the legend shows an N=4 entry.  Regenerate with distinct styling or note the coincidence in the caption.
- **A-4**: the MI values quoted to six decimals in supp.tex are only printed by `metric_compare.py`; the assertion window is ±2·10⁻³.  (My own recomputation confirms the printed digits.)
- **A-3**: `ns_general.py` Part A does not check p≥0 or S=2√2 for N=4..64 (my script does, for N=4..13 — passes).
- **A-5/A-6**: `review/chk_witness.py` has a check that can never fail and never exits non-zero.
- **A-7**: `make_fig_curves.py` transcribes the curve table by hand rather than computing it.
- **D-6**: `curves_highN.py` (exhaustive N=9..14, re-run clean) contradicts both papers' "N≥9 left to numerics, not performed" — the result exists and could be cited.
- **D-1/D-5/D-7/D-8/D-9**: PROOFS.md/RESULTS.md/REVIEW_REPORT.md/README use a theorem numbering that differs from the manuscript's (five theorems off, PROOFS mixes both); REVIEW_REPORT still lists as open several Tier-2 items and objections that were fixed; RESULTS §7 and both READMEs list 2 of 15 certificates and omit companion.tex.
- **A-2**: `verify.py` takes ~15 min, not "a few minutes".  **A-0**: `make_fig_curves.py` also writes an undeclared `paper/fig_curves.png`.

## 4. Minor / referee-facing suggestions (S3)

- **F5**: the non-monotone cap (1/4 at N=4, 2/5 at N=5) and the bends at F=1/N are artefacts of forcing a uniform baseline on all N states; with any baseline the cap is 1/4 for every N≥4.  Say so in one sentence, or a referee will.
- **F7**: the PR-box uniqueness lemma is standard — attribute (Popescu–Rohrlich 1994 is absent from refs.bib).
- **F12**: Mermin n=3 is combinatorially identical to CHSH; say so.
- **F11**: Theorem (Mermin n=5) could state equality F*=⌈3N/8⌉/N for residues r∈{0,1,2,3,4,6} mod 16, not only 0,1,2.
- **B2-f**: the H&B mutual-information formula in Table II equals theirs only for uniform settings (which they assume).
- **B12-d**: the Palmer rebuttal answers the coarse-graining argument but not his p-adic-metric argument.
- **B1-e**: make clear the "local determinism only" strengthening of Hall's bound is this paper's, not Hall's.
- supp.tex:26–28 "no floating point" over-claims (Parts B/C/G and the Garza–Hance script are float).
- 81 underfull hboxes; Table II lands a page late.

## 5. Agent findings I overrode

- Opus finding **C-3** ("`hp2020` defined but never cited") is a false positive: main.tex:94 cites it via `\onlinecite`, and Hossenfelder appears in the compiled bibliography.  Dropped.

## 6. Verification summary (what was actually run)

- `indep_check_chsh.py`: 100+ assertions, ALL PASS (runtime ≈ 4 min; LP-heavy).
- `indep_check_companion.py`: ALL PASS (runtime ≈ 14 min; Mermin-5 N=6 enumeration of 2.3 M multisets).
- `cglmp_quantum.py`: reproduces CGLMP 2.872934 and ADGL 2.914854.
- Opus agent 2 re-ran all 15 project certificates (logs in `rerun/`, `RUNSUMMARY.txt`): all exit 0, byte-identical to stored logs modulo timings; recompiled all three tex files into `build/` (exit 0; only the duplicate-label defect).
- Opus agent 1: 15 Crossref + 20 arXiv API/abs fetches; local PDFs read for Hall 2010, H&B 2020, Todorov, Garza–Hance, Palmer 2024, Rossi–Souza, CGLMP, Alai 2026b.  No fetch failed.

## 6b. Corrected copy of the papers (added later on 2026-09-03)

At the author's request the whole `paper/` directory was copied to `paper_fixed/` and every
S1/S2 fix above (plus the S3 clarifications) was applied there by `patch_paper.py`
(61 asserted exact-match edits); the original `paper/` is untouched.  Both papers now end
with an unnumbered "Note on provenance and corrections" explaining that the work was
produced by the local Qwen3.8-27B model under the human author's direction and was reviewed
and corrected by Claude Fable 5.1, listing the corrections.  `paper_fixed/CHANGELOG.md`
itemises every change; `paper_fixed/build/` holds the logs.  All three documents compile
with zero undefined references, zero duplicate labels and zero overfull boxes.  The
companion's CGLMP fix adds Acín–Durt–Gisin–Latorre's value as a second target; its N=2,3
values rest on `cglmp_Tprime.py` (floating-point exhaustive enumeration, labelled as such).

## 7. Process notes

- Pinned-model agent definitions were created at `~/.claude/agents/opus-high.md` and
  `~/.claude/agents/fable-high.md`, with usage notes in `~/.claude/CLAUDE.md`.  They are not
  loaded into an already-running session, so this review used `general-purpose` agents with an
  explicit Opus model override and effort instructions in the prompt; the pinned types will be
  available from the next session.
- Two Opus agents were used (sequentially, per the project's 5-slot LLM budget); no Fable agent
  was needed because the mathematical re-derivation was done in the main session.

## 8. Addendum (later on 2026-09-03): literature-positioning gap

While assessing the paper's contribution I noticed that neither paper cites Barrett & Gisin
(PRL 106, 100406, 2011), Pütz, Rosset, Barnea, Liang & Gisin (PRL 113, 190402, 2014), Pütz &
Gisin (NJP 18, 055006, 2016) or Koh et al. (PRL 109, 160404, 2012) — the central prior results on
how much measurement independence is needed, and the ones that establish that the answer is
strongly metric-dependent.  The project's two literature sweeps (Qwen's and my Opus agent's,
which only checked claims about *cited* works) both missed this.  All four references were
verified via Crossref and the arXiv API and added to `paper_fixed/` with a positioning
paragraph (`patch_paper2.py`, see `paper_fixed/CHANGELOG.md`).  Severity: S2 (a referee in
this subfield would certainly raise it).

## 9. Addendum: abstracts rewritten (user-approved) to lead with the contribution and scope; see paper_fixed/CHANGELOG.md third pass.

## 10. Addendum: residual points closed — N=9..14 curves independently re-derived (ALL PASS); T_Q confirmed as the I3 quantum maximum to numerical precision from ADGL 2002 + NPA 2008 (rank loop); companion wording strengthened; npa2008 added.

## 11. Addendum: supplement merged as Appendices A–D; placeholders filled (noaffiliation, acks, concurrent-submission cross-cites, code URL); cosmetics; standalone release folder bell-chsh-certificates/ built (see paper_fixed/CHANGELOG.md fifth pass). Items for the author to confirm: e-mail address, release repo name.
