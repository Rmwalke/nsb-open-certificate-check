# NS-B: an inspectable certificate companion

**Paper:** Global heat-seeded Picard convergence near a nonlinear periodic Navier–Stokes flow  
**Author:** Rushikesh Murlidhar Walke, Independent Researcher, Maharashtra, India  
**Version:** 1.1, 17 October 2026  
**Paper and published source:** https://doi.org/10.5281/zenodo.23429477

This companion makes one finite part of the preprint easier to inspect and rerun.

[Read the paper](Paper.pdf) · [Inspect the checker](audit/check_certificate.py) · [View the 17-block results](audit/results/W15_blocks.csv) · [Read the audit limits](audit/README.md) It is not an independent verification of the full proof, a peer-review report, or a novelty assessment.

## What the paper claims

For the unforced, viscosity-one Navier–Stokes equations on the normalized 2π-periodic torus, the manuscript considers smooth, real, mean-zero, divergence-free initial data W+h, where W=(sin z,sin x,sin y) and A_-1(h)≤10^-4. It claims convergence of ordinary full-data heat-seeded Picard iteration in X=L²((0,∞);A_0), with an error bound 6(9/10)^(n−3) for n≥3.

These are limited hypotheses. The claim does not address arbitrary smooth 3D data. The manuscript's separate enlarged regularity classes and supplementary fixed-W bounds are different claims.

## First question to inspect

In Proposition 2.4 / Section 2.5, do all 17 time-block inequalities for invariant-set preservation and third-iterate entry hold with the stored rational inputs?

The new script in `audit/` reads archived inputs as data, uses exact rational arithmetic, and does not import or execute the original scientific programs. Its result is conditional on the validity of the imported enclosures. See `audit/README.md` for the precise command, inputs and limitations, and `audit/results/` for the newly generated results.

The original `Code/` material and `Paper.pdf` were copied from the publication-v1.1 payload in the author's supplied handoff. Their bytes have not been edited. The live Zenodo record was not accessible to the retrieval tool used during preparation; no fresh byte comparison to the live deposit is claimed.

## What a successful finite check does not establish

It does not establish that every stored number bounds the intended Fourier field, signed residual, heat kernel or infinite tail. Those enclosure derivations require inspection. It also does not verify the analytic bilinear estimate, the full Galerkin-to-PDE argument, all manuscript claims, originality, or significance.

The article's Appendix A reports an exceptionally small positive invariance margin. Rounded decimal displays cannot decide its sign; use the exact fractions provided.

## Running and reviewing

The audit needs Python 3 and its standard library. It can run locally without a hosting account. `Inspect_Certificate.ipynb` is an optional notebook entry point, with a Binder environment specification. The notebook cells were executed locally. The hosted Binder build has not been tested. A proposed launch URL is https://mybinder.org/v2/gh/Rmwalke/nsb-open-certificate-check/main?urlpath=lab/tree/Inspect_Certificate.ipynb .

For the original, much broader reproduction instructions, read `Code/README.txt`. Replaying the original programs checks reproducibility and is distinct from this new, limited data-only arithmetic audit.

Useful feedback identifies the exact inequality, definition, enclosure, theorem or line being challenged and gives a derivation or counterexample. A result showing a real failure is useful and should be retained alongside the successful checks.

## Provenance and AI disclosure

The paper discloses substantive AI assistance in Section 5. This companion documentation, notebook and new audit code were also prepared with AI assistance on 18 October 2026. They have not received independent human review. Research authorship is as stated in the original paper.

Original license scopes and citation information are in `Code/LICENSE_SCOPE.txt` and `Code/CITATION.cff`. Preserve those files and notices when sharing the original material. The new companion files are supplied for inspection and reproducibility; they do not change the licensing or status of the original publication.
