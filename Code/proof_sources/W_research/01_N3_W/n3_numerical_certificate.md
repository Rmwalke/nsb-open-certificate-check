# N=3 certificate for the fixed datum W

2026-09-29. Private research record. This document describes a computer-assisted certificate, including the arithmetic trusted by its conclusion. It is not a formal proof-assistant certificate, external peer review, or a claim of novelty. The arbitrary-data three-dimensional Navier–Stokes problem remains open.

## Statement and conventions

Work on the torus `(R/(2*pi Z))^3` with normalized measure `(2*pi)^(-3) dx`, viscosity `nu=1`, zero spatial mean, and the Leray projection. Let

\[
 W=(\sin z,\sin x,\sin y),\qquad
 B(u,v)=\mathbb P((u\cdot\nabla)v).
\]

The Fourier cutoff is spherical: `P_N` retains the nonzero integer modes with `|k| <= N`. Let `v=v_3` solve the exact finite Galerkin system

\[
 v'+Av+P_3B(v,v)=0,\quad A=-\Delta,\quad v(0)=W.
\]

The residual is the **entire** exterior field `r=Q_3 B(v,v)`, with every generated output mode included. Under the validated enclosures recorded here, the usual H1 continuation theorem implies a unique global smooth solution with initial datum W. The same certificate implies global regularity for smooth, mean-zero, divergence-free data `W+h` with `||grad h||_2 <= 1/2500`; no Fourier-support restriction is imposed on h.

## Exact finite algebra

`n3_exact_algebra_intermediate.json` stores ten real divergence-free basis fields. A row with coefficient `a` at mode k means the complex Fourier coefficient `i*a`. The orthogonal basis has:

| Basis index | Laplacian eigenvalue | Squared L2 norm |
|---:|---:|---:|
| 0 | 1 | 3/2 |
| 1 | 2 | 3/4 |
| 2 | 5 | 3/16 |
| 3 | 5 | 15 |
| 4 | 6 | 1/8 |
| 5 | 6 | 12 |
| 6 | 8 | 6 |
| 7 | 9 | 1566 |
| 8 | 9 | 87/8 |
| 9 | 9 | 3/2 |

The first fields are W and `C=(sin y cos z, sin z cos x, sin x cos y)`. Exact rational convolution verifies that the span is invariant under the retained quadratic field. In coordinates,

\[
 x'_k=-\lambda_k x_k-\sum_{i\le j}c_{kij}x_ix_j,
 \qquad x(0)=(1,0,\ldots,0).
\]

For `i<j`, the stored coefficient already includes both ordered products; there is no further factor of two. There are 133 nonzero ODE coefficients. The 52 nonzero exterior pair fields generate 648 signed exterior modes, and their full Gram expansion gives 285 nonzero quartic terms in `||r||_2^2`. The union of retained basis supports has 108 signed modes and `sum |k|^2=660`. Closure is checked in the entire spherical cutoff, not only on those 108 modes.

The input SHA-256 is:

`15be3ffba859f976240e9c20a979d1a03a494dd69b8d404b873185bcfb9107fd`

## Finite-time enclosure

The primary run uses 1024 steps on `[0,8]`, step `h=1/128`, Taylor order 16, and outward integer intervals with scale `2^224`. An encoded pair `[a,b]` means `[a/2^224,b/2^224]`. All certificate decisions use integers and rational arithmetic. Floating-point numbers are only progress displays.

For each step, strict Picard inclusion

\[
 X+[0,h]F(\mathcal B)\subset\operatorname{int}\mathcal B
\]

first encloses the complete exact trajectory starting in X. The polynomial Taylor recurrence on X supplies initial derivatives; the same recurrence on the enclosing box supplies universal derivative remainders. Endpoints use variation of constants for each dissipative diagonal entry:

\[
 e^{-\lambda h}X+
 \sum_{n=0}^{p-1}w_{\lambda n}N_n(X)
 +w_{\lambda p}N_p(\mathcal B),\qquad
 w_{\lambda n}=\int_0^h e^{-\lambda(h-s)}s^n\,ds.
\]

An ordinary Taylor enclosure covers each entire time slab. Intersecting proven enclosures is allowed; replacing an enclosure by its midpoint is not done. Exact alternating series bound the weights and exponentials after range reduction; integer square roots give outward square-root bounds. A failed inclusion or an empty intersection rejects the run.

Let `G(t)=||grad v(t)||_(infinity,Frobenius)`. The first two basis fields obey

\[
 \|\nabla(aW+bC)\|_\infty=\sqrt{3(a^2+b^2)}.
\]

The remaining fields are bounded by the triangle inequality and their exact Fourier `sum |k| |hat b_k|`. Pointwise evaluation at the origin and the L2 gradient norm give valid lower bounds on G. The exterior residual norm is enclosed both by its exact quartic polynomial and by a midpoint value plus an L2 radius bound; their intersection remains an enclosure.

On every slab the code encloses and integrates

\[
 \Lambda(t)=\int_0^t(2G+2G^2)\,ds,\quad
 L(t)=\Lambda(t)-t/2,\quad
 \eta_T=\int_0^T2e^{-L}\|r\|_2^2\,dt,\quad
 I_T=\int_0^T e^{2L}\,dt.
\]

Simple rectangle bounds with full-slab intervals are sufficient; small endpoint uncertainty alone is not used as a quadrature error bound. Selected accumulated records and the final state are saved in `n3_validated_result.json`. Re-running the supplied code reconstructs the intervening steps.

## Analytic tail from time 8 to infinity

Let `S=8` and let E be a validated upper bound for `||v(S)||_2^2`. Exact Galerkin energy and Poincare give

\[
 \|v(S+\tau)\|_2^2\le Ee^{-2\tau}.
\]

Cauchy–Schwarz on the retained support and the cutoff bound give

\[
 G(S+\tau)\le H e^{-\tau},\qquad H=\sqrt{660E},
\]

\[
 \|r(S+\tau)\|_2^2
 \le\|v\|_\infty^2\|\nabla v\|_2^2
 \le108\cdot9 E^2 e^{-4\tau}=972E^2e^{-4\tau}.
\]

Consequently `Lambda(infinity)-Lambda(S) <= Delta=2H+H^2`, and

\[
 \eta_{\rm tail}\le\frac47\,972E^2e^{-L(S)},\qquad
 I_{\rm tail}\le e^{2L(S)+2\Delta}.
\]

These are upper bounds only. The residual projection cannot increase its L2 norm. The finite Galerkin flow exists throughout this tail by its energy identity, so the tail does not assume existence of the unknown infinite-dimensional solution.

## Recorded result and strict rational margins

The primary run completed `[0,8]` with status `PASS_ORIGINAL_GATE`. It gives

\[
 \kappa=27(25/3)^6\eta^2 I
 \le
 \frac{13108650706211112106568919404571274711326726197894008841903862847223}
 {26959946667150639794667015087019630673637144422540572481103610249216}
 <0.49.
\]

The decimal value of the displayed upper bound is approximately 0.486226878267. Convenient outward rational bounds are

\[
 \eta<\frac{587}{10^9},\qquad I<156400,\qquad
 \Lambda(\infty)<\frac{176}{25},\qquad e^{176/25}<1150.
\]

Even substituting these rounded bounds gives

\[
 27(25/3)^6(587/10^9)^2\,156400<49/100<1.
\]

The primary gate passes without needing the alternative forcing-time criterion also recorded by the runner. `n3_pde_scalar_audit.md` supplies the H1 error inequality, scalar comparison, all-time continuation, and the consequences

\[
 \|\nabla(u_M-v)\|_2^2<10^{-3}e^{-t/2}\quad(M\ge3),
 \qquad
 \sup_M\sup_{t\ge0}\|\nabla u_M(t)\|_2^2<14,
\]

\[
 \sup_M\int_0^\infty\|\nabla u_M(t)\|_2^4\,dt<21/2.
\]

The same statements for W concern exact Galerkin solutions and their PDE limit, not the numerical center trajectory. The perturbation-ball conclusion is separate; these W-specific energy bounds are not assigned unchanged to perturbed data.

## Sobolev input, with normalization checked

The proof uses the vector/tensor estimate `||f||_6 <= (25/3)||grad f||_2` for componentwise mean-zero f. Here is a short derivation, rewritten from the attached written Sobolev note without executing its historical code.

For any nonnegative smooth scalar a on the period-L torus, variation on the two arcs between a coordinate-circle minimum and maximum gives a bound by its coordinate mean plus one half of the unnormalized derivative integral. Multiplying the three coordinate bounds and applying repeated Cauchy–Schwarz yields

\[
 \|a\|_{3/2}\le\|a\|_1+\frac L6\sum_i\|\partial_i a\|_1.
\]

Apply this to `a=|f|^4`, use Holder, and divide by `||f||_6^3` when f is nonzero:

\[
 \|f\|_6\le\|f\|_2+\frac{2L}{\sqrt3}\|\nabla f\|_2.
\]

At `L=2*pi`, mean-zero Poincare and `pi<22/7`, `sqrt(3)>12/7` give `1+4*pi/sqrt(3)<25/3`. The zero field is immediate. The same argument applies to the gradient tensor, which is componentwise mean zero; Parseval gives `||grad grad f||_2=||Delta f||_2`. No target-dimension factor is needed.

## Verification scope

The exact Fourier algebra was independently reconstructed, and the numerical core and diagnostics were read separately. A second finite run uses ordinary Taylor endpoints at order 32 and 384-bit outward intervals; it completed `[0,8]` and also passed the original gate with an outward upper bound of `0.486226878267`. All 330 coordinate enclosures at the 33 common recorded times intersect. That run shares the interval arithmetic and diagnostics, so it is an alternative propagation check, not a wholly separate software implementation. A separate rational positive-series calculation also reconstructs the primary analytic-tail bound.

The certificate still trusts the supplied Python implementation, its execution environment, the written analytic estimates, and standard local existence/H1 continuation. A source hash identifies bytes; it does not prove mathematics. No historical script, test, frozen verifier, or import was executed for this work. No conclusion is asserted for arbitrary initial data, arbitrary viscosity, or a differently normalized torus.
