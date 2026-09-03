# Superseded outputs of the original search pipeline

These files record how the results were *first found* during the Qwen phase (floating-point LP
bisections, feasibility probes, an embedding check, a self-review pass, the figure script's stdout).
All of them were later replaced by the exact-arithmetic certificates in `../` and are kept only as
history.  In particular the bisection values `ns_bisect_N{3,4}.out` (F* ≈ 0.13812, ≈ 0.10358) are
resolution artefacts of a 13-step bisection, as the main paper's Appendix C explains; the exact values
are (√2−1)/3 and (√2−1)/4.  `referee_pass.out` is a self-review transcript written by the Qwen model
before the independent review.  Nothing in the papers depends on these files.
