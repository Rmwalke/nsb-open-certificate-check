# NS-B: an inspectable certificate companion

> **Reviewer quick start — a limited, conditional arithmetic audit**
>
> The [recorded Python audit](audit/results/W15_data_only_audit.json) reports `PASS_DATA_ONLY`: **17/17** rectangle-invariance comparisons, **17/17** strict third-iterate entry comparisons, and **51/51** archived prefix recurrence comparisons. These are stored results; the documentation review on 28 October 2026 did not rerun the research calculations.
>
> Acceptance uses Python `Fraction`, not rounded decimals. The smallest recorded invariance margin is approximately **1.1507455879810338 × 10⁻¹⁹**; its exact fraction is in `summary.smallest_invariance_margin_exact`. See [the checker](audit/check_certificate.py) and [its scope](audit/README.md).
>
> The checker reads **three** archived JSON inputs: the W15 perturbation gate and prefix result, plus the W14 refined baseline. It checks their dependency hashes and consistency, then recomputes the finite inequalities. Those internal hash checks do not independently authenticate a release or prove the imported analytic enclosures.
>
> A successful check is conditional on those enclosures. It does **not** verify the complete PDE proof, establish novelty, or constitute peer review. The manuscript's claim does not address arbitrary smooth 3D data. The original `Code/` bundle and `Paper.pdf` are preserved unchanged.

**Check in your browser:** [NS-B Experiment Lab](https://nsb-experiment-lab-rmwalke.rushikeshwalkehtc.chatgpt.site) recomputes **34 selected comparisons** from a derived input snapshot, using BigInt fractions for decisions. It does not perform the Python audit's additional archived recurrence, original-file hash and consistency checks; browser PASS is not the complete `PASS_DATA_ONLY` result. Its arithmetic source was inspected in this review; no new browser calculation was run.

<details>
<summary>Exact audit inputs and license scope</summary>

All paths are relative to this repository:

- `Code/proof_sources/W_research/15_perturbed_picard/perturbation_block_gate_result.json`
- `Code/proof_sources/W_research/15_perturbed_picard/perturbed_picard_prefix_static_result.json`
- `Code/proof_sources/W_research/14_unit_W_picard/refined_W_block_gate_result.json`

The root [MIT license](LICENSE) applies to the new companion software in `audit/check_certificate.py` and the code cells of `Inspect_Certificate.ipynb`. It does not relicense `Paper.pdf`, certificate data or the original publication bundle. The original distribution retains the MIT/CC BY 4.0 scopes specified in [Code/LICENSE_SCOPE.txt](Code/LICENSE_SCOPE.txt), including its attribution notices.

</details>

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
