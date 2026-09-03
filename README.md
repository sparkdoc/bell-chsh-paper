# bell-chsh-paper

## What this repository is: a physics paper written by a small local language model

This repository is the record of an experiment: **can a small, locally hosted open-weight language
model do research-level theoretical physics, from choosing a problem to a submittable manuscript?**
The model was asked to identify a scientific problem of interest, research the literature, do the
mathematics and computations, and write a publication.  The result is two papers on the minimum
measurement dependence ("superdeterminism cost") needed to reproduce quantum Bell-inequality
violations:

* **Exact minimum measurement dependence for violating the CHSH inequality** (`paper/main.tex`,
  with the former supplement folded in as Appendices A–D), and
* **Fine-tuning beyond CHSH: minimum measurement dependence for CGLMP, Mermin, and contextuality
  inequalities** (`paper/companion.tex`).

### The papers

| | Qwen3.8-27B original (unreviewed, 2026-09-02) | after review and correction by Claude Fable 5.1 (2026-09-03) |
|---|---|---|
| **Main paper** — *Exact minimum measurement dependence for violating the CHSH inequality* | [PDF](paper/original_qwen_2026-09-02/main.pdf) (12 pp.) + [Supplemental Material PDF](paper/original_qwen_2026-09-02/supp.pdf) (3 pp.) | [PDF](paper/main.pdf) (16 pp., supplement folded in as Appendices A–D) · [LaTeX](paper/main.tex) |
| **Companion paper** — *Fine-tuning beyond CHSH: minimum measurement dependence for CGLMP, Mermin, and contextuality inequalities* | [PDF](paper/original_qwen_2026-09-02/companion.pdf) (10 pp.) | [PDF](paper/companion.pdf) (10 pp.) · [LaTeX](paper/companion.tex) |
| **What changed between the two columns** | | [`paper/CHANGELOG.md`](paper/CHANGELOG.md) (reader-facing summary) · [full review report](review_2026-09-03/REVIEW_REPORT_fable_2026-09-03.md) |

The left column is the unaltered output of the local model; the right column is what will be
submitted to arXiv.  The mathematics is the same in both; the differences are a wrong quantum
target in the companion, bibliography errors, missing prior literature, and wording — see the
changelog.

The human author directed the process (choice of direction, approval of actions, session
management) and is the author of record, but supplied minimal scientific or technical input.  On 2026-09-03 the finished manuscripts and the whole verification pipeline were critically
reviewed by a second, much larger model, Claude Fable 5.1, which confirmed the mathematics and
corrected the defects it found (bibliography, one wrong cross-reference, the CGLMP quantum target,
literature positioning, wording).  The papers as they stood **before** that review are preserved
unaltered in `paper/original_qwen_2026-09-02/`; the corrected papers are in `paper/`; every change is
itemised in `paper/CHANGELOG.md` and in `review_2026-09-03/REVIEW_REPORT_fable_2026-09-03.md`.  Each
paper ends with a "Note on provenance and corrections" saying the same thing in print.

### The setup

* **Model.** Qwen3.8-27B in the Unsloth GGUF release, quantisation `UD-Q8_K_XL`
  (`unsloth/Qwen3.8-27B-GGUF:UD-Q8_K_XL`), served locally by Unsloth Studio through an
  OpenAI-compatible endpoint with five concurrent generation slots.
* **Hardware.** Two NVIDIA RTX 6000 Ada GPUs (48 GB each).
* **Harness.** Anthropic's Claude Code agentic CLI, pointed at the local endpoint instead of a
  hosted model.  Claude Code supplied the tool loop (shell, file editing, web fetching), sub-agents,
  multi-agent workflows, persistent memory and the permission system; the local model supplied every
  token of reasoning and text.  The five slots were shared by the main session, sub-agents and the
  harness's own auxiliary calls (e.g. the command-safety classifier), which shaped how much
  parallelism was possible (see `CLAUDE.md`-style operating notes quoted in the review report).
* **Verification discipline.** Every exact claim was required to carry a re-runnable certificate
  (exact rational or quadratic-field arithmetic, exit code 0 iff all assertions pass); those
  certificates are in `certificates/` with their logs in `logs/`.
* **Review.** The same harness, driving Claude Fable 5.1 with Claude Opus 5 sub-agents, performed
  the independent review; its from-scratch checks are in `review_2026-09-03/`.

### How much generation it took

Token counts below come from the harness's stored transcripts: the main interactive sessions and
the 178 sub-agent and multi-agent-workflow transcripts that the three main Qwen working sessions
persisted alongside them (`review_2026-09-03/analyze_transcripts.py <project-slug>` → `out_transcripts.txt`;
usage de-duplicated per message id).  A few small early sessions kept no sub-agent transcripts, so
the Qwen figures are still lower bounds, but close ones.  "Prompt tokens" is the total context the
model re-read across all turns (uncached input plus cache reads); the local server reported no
prompt caching for sub-agents.

| phase | dates | assistant messages | generated (output) tokens | prompt tokens processed | tool calls |
|---|---|---|---|---|---|
| Qwen3.8-27B, main sessions | 2026-08-21 → 2026-09-03 | ≈ 1,540 | ≈ 2.24 M | ≈ 187 M | ≈ 2,040 |
| Qwen3.8-27B, sub-agents and workflows (178 transcripts) | 2026-08-21 → 2026-09-03 | ≈ 2,850 | ≈ 6.59 M | ≈ 275 M | ≈ 3,340 |
| **Qwen3.8-27B, total** | 13 calendar days | **≈ 4,390** | **≈ 8.8 M** | **≈ 460 M** | **≈ 5,400** |
| Claude Fable 5.1 + 2 Opus 5 sub-agents (review and corrections) | 2026-09-03 | ≈ 260 | ≈ 0.36 M | ≈ 52 M (almost all cache reads) | ≈ 280 |

Most of the Qwen sub-agent volume came from one session: the adversarial deep review of
2026-08-30, a multi-agent workflow plus follow-up re-derivations, accounts for 87 of the 178
transcripts and about 5.5 M generated tokens.

**Decode time.** On this hardware the server's aggregate decode rate was roughly 65 tokens/s,
shared among whatever agents were generating concurrently (concurrency did not raise the total).
The 8.8 M generated tokens therefore correspond to about 38 hours of pure generation, a floor
that no amount of parallelism could lower.

**Prefill (KV-cache load) time — an estimate.** Summed over every turn, the model's context
totalled ≈ 462 M tokens, but with prefix (KV-cache) reuse only the new material added each turn
has to be processed; the recorded context lengths put that incremental volume at ≈ 18 M tokens
(4 % of the total).  At the ≈ 1,000–2,000 tokens/s prefill typical of a 27B 8-bit model on this
class of GPU, that is about 2.5–5 hours.

**Bottom line.** Under normal operating conditions (KV caching on), producing these two papers with
Qwen3.8-27B on two RTX 6000 Ada GPUs takes **roughly 40–45 hours of GPU time**: ≈ 38 h of decoding
at 65 tokens/s plus 2.5–5 h of prefill, spread here over 13 calendar days of intermittent sessions.
To estimate a similar project on your own hardware, divide the ≈ 8.8 M generated tokens by your
decode rate (e.g. at 20 tokens/s, about 120 hours) and the ≈ 18 M incremental prompt tokens by
your prefill rate.

### What the experiment showed

The model's mathematics held up: every theorem, curve, witness and certificate the reviewer could
test was correct.  Its weaknesses were bibliographic (four wrong or fabricated reference details,
two of them annotated in-file as "verified"), a missed literature thread (Barrett–Gisin 2011, Pütz
et al. 2014), one wrong quantum target (the CGLMP maximally-entangled value taken for the maximum),
and a number of wording and record-keeping inconsistencies.  The most telling finding came from
following up the missed thread: the paper's headline for CHSH — S ≤ 2+8F with threshold
F* = (√2−1)/4 — turned out to be a re-parametrisation (ℓ = 1/4 − F) of the measurement-dependent-
locality bound 4(1−2ℓ) that Pütz et al. published in 2014, built on the same extremal model.  The
model derived a known result by a different route, in a new variable, without recognising it.
Whether that reflects latent recall of training data or independent rediscovery cannot be
determined from the outside; either way, it is the failure mode to expect from LLM-generated
research, and the reason the novelty claims are stated as search results rather than as facts.  Scientifically the papers are a
correct, modest, incremental contribution: an exact, certified answer in one natural metric to a
question whose qualitative answer was already known.  The full assessment is in the review report.

## Layout

| directory | contents |
|---|---|
| `certificates/` | the exact-arithmetic certificates and LP cross-checks cited in the papers (Python 3) |
| `logs/` | the certified output of every script, as cited (numbers quoted in the papers are taken from these); `logs/lit2/` holds the raw arXiv-API responses and fetch scripts behind the literature sweep, `logs/superseded/` the outputs of the original search pipeline that the exact certificates replaced |
| `witnesses/` | saved witness arrays `solN{2,3,4}_{p,u,v}.npy`, `ns_solN{3,4}_{p,u,v}.npy`, `lg_results.json` |
| `review_2026-09-03/` | the independent review's from-scratch verification scripts and outputs, and the review report |
| `paper/` | LaTeX sources, bibliography, figure and compiled PDFs of both papers; `paper/original_qwen_2026-09-02/` holds the unaltered pre-review PDFs |

## Running

Requirements: Python 3.8+; `numpy` and `scipy` (HiGHS) only for the floating-point LP cross-checks;
`matplotlib` only for `make_fig_curves.py`.  Every script resolves its paths relative to this
repository, so it can be run from any working directory:

```
python3 certificates/audit1.py          # exact S_max curves N=2..8, N=4 witness in Q(sqrt2), NS N=2  (~1 min)
python3 certificates/curves_highN.py    # exact curves N=9..14                                         (~30 s)
python3 certificates/ns_general.py      # NS witnesses, all N>=3
python3 certificates/prbox_ns.py        # NS algebraic maximum (PR-box decomposition)
python3 certificates/onesided.py        # one-sided measurement dependence
python3 certificates/metric_compare.py  # metric conversions, mutual information
python3 certificates/fringe_pm.py       # Peres-Mermin square (+ aggregation robustness)
python3 certificates/ineq_cglmp.py  certificates/ineq_mermin.py  certificates/ineq_mermin5_n5.py  certificates/ineq_xcheck.py
python3 review_2026-09-03/indep_check_chsh.py         # independent re-verification of the main paper   (~4 min)
python3 review_2026-09-03/indep_check_companion.py    # independent re-verification of the companion    (~15 min)
```

Each exact certificate exits with code 0 if and only if every assertion passes.  Scripts write their
logs into `logs/`, overwriting the shipped copies; the shipped copies are the record cited in the
papers (`verify.py` takes about 15 minutes; `phase2.py`, `phase3.py`, `capregion.py`, `ns_tight.py`
are the original LP search pipeline and are slow).

### Which script certifies what

| claim (main paper) | script → log |
|---|---|
| exact curves S_max(N,F), N=2..8; cap ⌈N/4⌉/N; N=4 witness in Q(√2); NS N=2 flat | `audit1.py` → `audit1.out` |
| exact curves N=9..14 | `curves_highN.py` → `curves_highN.out`; independently `review_2026-09-03/curves_N9to14.py` |
| NS witnesses N≥3 (Thm 7) | `ns_general.py` → `ns_general.out`; `ns_witness.py`/`ns_verify.py` → `ns_witness.out` |
| NS algebraic maximum (Thm 8) | `prbox_ns.py` → `prbox_ns.out` |
| one-sided measurement dependence (Prop. 2) | `onesided.py` → `onesided.out` |
| metric comparison, MI values (Sec. VII.5) | `metric_compare.py` → `metric_compare.out` |
| aggregation robustness (Prop. 3) | `fringe_pm.py` Part PM-7 → `fringe_pm.out` |
| Garza–Hance figures (Sec. VII.4) | `chk2_superdet_garza.py` → `out_superdet_garza.txt` |
| LP cross-checks | `verify.py` → `verify_run.log`; `capregion.py` → `capregion.out`; `ns_tight.py` |
| **companion**: CGLMP / Mermin / Peres–Mermin | `ineq_cglmp.py`, `ineq_mermin.py`, `ineq_mermin5_n5.py`, `ineq_xcheck.py`, `fringe_pm.py` → same-named `.out` (`ineq_xcheck2.out`) |
| CGLMP at the Acín et al. target | `review_2026-09-03/cglmp_quantum.py`, `cglmp_Tprime.py` |

## Provenance of the code

`certificates/*.py` were written by the Qwen3.8-27B model during the project (2026-08); the only
change made for this release is that hard-coded absolute paths were replaced by paths relative to the
repository root (one header line per file).  `review_2026-09-03/*.py` were written from scratch by
Claude Fable 5.1 and share no code with the certificates.

## License

MIT (see `LICENSE`).  The license covers the code, logs and documentation in this repository; the
manuscripts in `paper/` remain the author's copyright pending publication.
