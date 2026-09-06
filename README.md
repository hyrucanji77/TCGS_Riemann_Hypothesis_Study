# Moment and Prime-Trace Criteria for the Riemann Hypothesis, with Two Spectral Exclusions

**Henry Arellano-Peña** · NEBIOT S.A.S., Bogotá, Colombia  
**Manuscript version:** SUBMISSION_FINAL · **Manuscript date:** 6 September 2026

[Read the paper (PDF)](main.pdf) · [Complete LaTeX source](main.tex) · [Download the Overleaf bundle](TCGS_Riemann_Hypothesis_Study_Overleaf.zip)

## Scope

This repository contains the final submission manuscript and its reproducible computational diagnostics. The paper develops moment and prime-trace criteria, finite-to-infinite positivity transfer, and spectral exclusions, including the constant-coefficient 2:1 Morse family with the stated real constant separated endpoint conditions and constant spectral shifts.

**The Riemann hypothesis is not proved.** The source-derived arithmetic realization is identified as an open construction. The TCGS–SEQUENTION realization problem is distinguished from the unconditional analytic results. Numerical reproduction is not an interval certificate, a proof of RH, or independent peer review.

## Paper and source integrity

`main.tex` is the unchanged final submission LaTeX file. Its SHA-256 is:

```text
ac81d0f3aaa5cbe1df54108a6165efcee2e29885637b9c21aa32d6dbbb490df6
```

The repository PDF is compiled on GitHub from that source. It is not represented as a byte-identical upload of an earlier PDF: PDF metadata and the TeX environment can change the binary file. Build checksums, TeX version, page information, and extracted paper text are recorded under `verification/`.

The repository preserves the original nine diagnostic scripts, the release-verification driver, the pinned requirement, and all thirteen reference JSON outputs. Repository documentation and build automation are additional delivery files; they do not amend the manuscript. This repository is the paper-and-code release, not the separate full historical/source-PDF archive. It does not include private drafting correspondence, superseded manuscripts, or the generated infographic.

## Compile the paper

Upload the Overleaf bundle and select the root-level `main.tex`. Use **pdfLaTeX**. The numbered bibliography is embedded; no separate `.bib` file or Python execution is needed to compile.

For a local TeX installation:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The GitHub workflow in `.github/workflows/build-paper.yml` builds `main.pdf` and `TCGS_Riemann_Hypothesis_Study_Overleaf.zip`. The ZIP is regenerated from the tracked source, PDF, code, and reference outputs rather than copied from a historical release archive.

## Reproduce the computations

The original reference outputs were generated with **Python 3.13.5** and **mpmath 1.3.0**. A matching Python version permits byte-for-byte comparisons including environment metadata. The verifier separately reports equality with Python-version fields ignored.

```bash
python -m pip install -r requirements.txt
python verify_release.py --output-dir verification/fresh_run
```

Choose a new output directory for every run. The verifier copies only the nine scripts and `requirements.txt` into its fresh working directory. It then runs thirteen diagnostic configurations, compares all generated JSON files with the reference outputs, and checks sixty numerical values against Tables 1–2 and Equations (13.20)–(13.21) in `main.tex`. It writes logs and `report.json` in the requested directory. Expected JSON files are not used as computational inputs.

| Script | Purpose |
|---|---|
| `diagnostics.py` | Central xi moments, Hankel determinants, exponential comparison |
| `morse_diagnostics.py` | Calibrated Dirichlet Morse asymptotics |
| `robin_shift_diagnostics.py` | Robin normalization, shift asymptotics, exact coefficients |
| `boundary_trace_diagnostics.py` | Exact Schur-complement and noncommuting resolvent checks |
| `scale_endpoint_diagnostics.py` | Free scale/endpoint comparison and quarter-shift normalization |
| `shifted_tower_diagnostics.py` | Shifted Stieltjes prefix and even/odd block comparisons |
| `completed_tower_diagnostics.py` | Both moment towers through degree four using Taylor order twenty |
| `kinetic_rate_diagnostics.py` | Exact kinetic–rate reduction and boundary transport |
| `sesquilinear_diagnostics.py` | Exact Gaussian-rational conjugate-linearity checks |

`verify_release.py` is the common reproduction driver. Its generated report distinguishes runtime failures, file equality, printed-value agreement, and exact finite algebraic checks. Passing these tests verifies the declared computations, not the missing infinite-dimensional arithmetic identity.

## Citation and licensing

Citation metadata are provided in [CITATION.cff](CITATION.cff). No journal acceptance or repository DOI is claimed.

The manuscript and original code are licensed under **Creative Commons Attribution 4.0 International (CC BY 4.0)**, consistent with the manuscript. See [LICENSE](LICENSE). Third-party works retain their own licenses.

The computational-assistance disclosure in the manuscript is retained unchanged. Responsibility for mathematical verification and publication rests with the named author.
