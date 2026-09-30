# Moment and Prime-Trace Criteria for the Riemann Hypothesis, with Spectral Exclusions and Nonlocal Positivity

**Henry Arellano-Peña** · NEBIOT S.A.S., Bogotá, Colombia  
**Release:** REVIEW_AMENDED_2026-09-30 · **Original manuscript date:** 6 September 2026 · **Amended:** 30 September 2026

[Read the paper (PDF)](riemann_hypothesis_TCGS_SEQUENTION.pdf) · [Complete LaTeX source](main.tex) · [Download the Overleaf bundle](TCGS_Riemann_Hypothesis_Study_Overleaf.zip)

## Scope

This repository contains the reviewed manuscript and its reproducible computational diagnostics. The paper develops reciprocal-square moment and completed prime-trace criteria, finite-to-infinite positivity transfer, compact-elliptic and regular boundary-power exclusions, and the full constant-coefficient 2:1 Morse-family exclusion with the stated separated boundary conditions and constant spectral shifts.

Appendix E adds Weil positivity and screw kernels, the exact Laplace and characteristic-ratio normalizations, multiplicity bookkeeping, conditional source-to-kernel Gram transfer, nonlocal real-zero limits, and finite-window arithmetic error estimates. The title, abstract, versioned citations, and script inventory follow the review-amended manuscript.

**The Riemann hypothesis is not proved.** The source-derived arithmetic realization and the required nonlocal analytic limit remain unconstructed. The TCGS–SEQUENTION realization problem is distinguished from the unconditional analytic results. Numerical reproduction is not an interval certificate, an RH proof, or independent peer review.

## Paper and source integrity

The root-level `main.tex` is byte-identical to the approved review-amended source. Its SHA-256 is:

```text
5a6dc2b0ddc6752466a3d1b9b21411007f95cd98b88c29264cd3baaaf9456069
```

The repository's `main.pdf` and `riemann_hypothesis_TCGS_SEQUENTION.pdf` are identical aliases of the 65-page paper compiled on GitHub from that source. The rebuilt PDF is not claimed byte-identical to the local review PDF: PDF metadata and the TeX environment can change its binary checksum. The checksum of the approved local PDF is retained separately in `release_integrity_20260930.json`.

Build checksums, TeX version, page information, extracted text, the exact-source import record, and both numerical verification reports are recorded under `verification/github_*`. The final source and each scientific file are checked against the supplied release checksums before packaging. Repository history is retained.

This is the paper-and-code release, not the separate historical/source-PDF archive. It does not publish private correspondence, superseded manuscript copies, the third-party source-PDF collection, or the generated infographic. The complete active manuscript and its embedded bibliography are included. Original scientific filenames, the nine legacy scripts, the original verifier, and their thirteen reference outputs are preserved; the Appendix E script, its verifier, and two reference outputs are added.

## Compile the paper

Upload the Overleaf bundle, select the root-level `main.tex`, and use **pdfLaTeX**. The numbered bibliography is embedded; no separate `.bib` file or Python execution is needed for compilation.

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

`.github/workflows/build-paper.yml` compiles both PDF aliases, runs both verification drivers, and regenerates `TCGS_Riemann_Hypothesis_Study_Overleaf.zip`. The metadata-packaging workflow checks both successful reports before refreshing the bundle. `tools/publish_review.py` is repository delivery bookkeeping, not an additional scientific diagnostic or proof checker.

## Reproduce the computations

Reference outputs use **Python 3.13.5** and **mpmath 1.3.0**. The verifiers distinguish byte equality from equality after ignoring Python-version metadata.

```bash
python -m pip install -r requirements.txt
python verify_release.py --output-dir verification/fresh_legacy_run
python verify_screw_release.py --output-dir verification/fresh_screw_run
```

Use new, nonexistent output directories. There are **ten scientific diagnostic scripts and two verification drivers**. `verify_release.py` reruns nine scripts in thirteen configurations, checks sixty printed values, and covers 122 distinct exact finite checks. `verify_screw_release.py` separately runs the Appendix E script at 40 and 60 digits, checks the three rows of Table 5 at both precisions, and covers 56 distinct exact checks counted once across the precisions. Together they check fifteen reference outputs, 66 printed-value comparisons, and 178 distinct exact checks. Expected JSON files are compared only after fresh computation; they are not scientific input to the diagnostic calculations.

| Scientific script | Purpose |
|---|---|
| `diagnostics.py` | Central xi moments, Hankel determinants, exponential comparison |
| `morse_diagnostics.py` | Calibrated Dirichlet Morse asymptotics |
| `robin_shift_diagnostics.py` | Robin normalization, shift asymptotics, exact coefficients |
| `boundary_trace_diagnostics.py` | Exact Schur-complement and noncommuting resolvent checks |
| `scale_endpoint_diagnostics.py` | Free scale/endpoint comparison and quarter-shift normalization |
| `shifted_tower_diagnostics.py` | Shifted Stieltjes prefix and even/odd block comparisons |
| `completed_tower_diagnostics.py` | Both moment towers through degree four using Taylor order twenty |
| `kinetic_rate_diagnostics.py` | Exact kinetic–rate reduction and boundary transport |
| `sesquilinear_diagnostics.py` | 35 exact Gaussian-rational conjugation checks, including negative controls |
| `screw_kernel_diagnostics.py` | Appendix E finite arithmetic kernels, ratio normalizations, and exact error/counting checks |

Passing these tests verifies the declared finite computations, not the missing infinite-dimensional arithmetic identification. GitHub records the legacy and Appendix E results separately in `verification/github_diagnostics_report.json` and `verification/github_screw_diagnostics_report.json`.

## Citation and licensing

Citation metadata are provided in [CITATION.cff](CITATION.cff). No journal acceptance or repository DOI is claimed.

The manuscript and original code are licensed under **Creative Commons Attribution 4.0 International (CC BY 4.0)**, consistent with the manuscript. See [LICENSE](LICENSE). Third-party works retain their own licenses.

The manuscript's computational-assistance disclosure is retained unchanged. Responsibility for mathematical verification and publication rests with the named author.
