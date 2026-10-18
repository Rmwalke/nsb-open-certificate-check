# True heat-seeded Picard convergence near unit W

Private local written proof, 2026-09-30 (Asia/Kolkata). The new checker
uses earlier result JSON files only as data. It does not execute or
import earlier scripts. All time integrals below are analytic endpoint
moments; no ODE/PDE evolution, quadrature, web access, or upload is used.

**For smooth real mean-zero divergence-free initial data
\(u_0=W+h\), where \(W=(\sin z,\sin x,\sin y)\) and
\(A_{-1}(h)\le10^{-4}\), the ordinary Picard sequence starting from
the full heat flow \(e^{-tA}(W+h)\) converges on \([0,\infty)\)
in \(\mathcal X=L^2(0,\infty;A_0)\).** The finite-prefix check
here proves that this actual sequence enters the independently
certified invariant block ball after its third iterate. The datum
\(h\) may have arbitrary Fourier support. The resulting explicit
error bound is
\(\|v_h^{[n]}-u_h\|_{\mathcal X}\le6(9/10)^{n-3}\) for \(n\ge3\).

The separate regularity result already covered the larger class
\(A_{-1}(h)\le1/20\) at unit \(W\). The new statement concerns
convergence of the ordinary iteration on a smaller explicit
neighborhood. It is not an enlarged global regularity class or an
all-data convergence result.

## 1. The iteration and the retained initial perturbation

Use the normalized periodic viscosity-one setting,
\(A=-\Delta\), and
\(A_s(f)=\sum_{k\ne0}|k|^s|\widehat f_k|\). Let
\(\mathcal Q\) be the symmetric Navier--Stokes heat bilinear map,
with the established constant

\[
 \|\mathcal Q(f,g)\|_{\mathcal X}
 \le c\|f\|_{\mathcal X}\|g\|_{\mathcal X},
 \qquad c=\frac1{2\sqrt2}.
\]

Define the ordinary maps and their true heat-seeded sequences by

\[
 T_h(v)=e^{-tA}(W+h)-\mathcal Q(v,v),\qquad
 v_h^{[0]}=e^{-tA}(W+h),\qquad
 v_h^{[n+1]}=T_h(v_h^{[n]}).
\]

Write \(w_n=v_0^{[n]}\). Thus \(w_0=e^{-t}W\), and
\(w_1,w_2,w_3\) are exactly the previously constructed corrected
references. Every perturbed iterate uses the entire initial datum
\(W+h\). In particular, the perturbation is present both in the
seed and in every map application; it is not a zero-initial-value
residual subtraction from an unperturbed reference.

The 17 intervals \(I_i=[a_i,b_i]\) have endpoints

\[
 0,\tfrac1{16},\tfrac18,\tfrac3{16},\tfrac14,\tfrac5{16},
 \tfrac38,\tfrac7{16},\tfrac12,\tfrac58,\tfrac34,\tfrac78,
 1,\tfrac54,\tfrac32,2,3,\infty.
\]

Let \(\kappa_{ij}\) denote the rigorously enclosed causal heat
memory factors from the block certificate, with
\(\kappa_{ii}=1\). For \(j<i\), they satisfy

\[
 \kappa_{ij}\ge\sup_{\lambda\ge1}
 e^{-\lambda(a_i-b_j)}
 \sqrt{1-e^{-2\lambda(b_i-a_i)}}.                 \tag{1}
\]

On the infinite final block, the square-root factor is one. Source
blocks contribute at their original times through (1); no boundary
value or error history is reset.

## 2. Uniform block caps for the full heat perturbation

For \(\delta=A_{-1}(h)\), Fourier-mode Minkowski gives

\[
 \|e^{-tA}h\|_{\mathcal X(I_i)}
 \le\frac\delta{\sqrt2}\sup_{\lambda\ge1}
 e^{-\lambda a_i}\sqrt{1-e^{-2\lambda(b_i-a_i)}}
 \le s_i.                                         \tag{2}
\]

The perturbation block certificate encloses these factors by positive
rational numbers, evaluated at \(\delta=1/10000\). This uniformly
covers all smaller \(\delta\). In particular, on the first block
\(a_1=0\), the supremum is one: the correct uniform cap is
\(\delta/\sqrt2\), regardless of how short that first block is.
The proof does not assume an unjustified small-time gain for
arbitrarily high Fourier modes of \(h\).

Let \(V_{n,i}\ge\|w_n\|_{\mathcal X(I_i)}\), for \(n=0,1,2\).
The first two caps follow from

\[
 A_0(w_0(t))=3e^{-t},\qquad
 A_0(w_1(t))=3e^{-t}+3te^{-2t}.
\]

The corresponding squared tail moments are exactly

\[
 F_0(t)=\frac92e^{-2t},
\]
\[
 F_1(t)=\frac92e^{-2t}+(6t+2)e^{-3t}
 +\left(\frac94t^2+\frac98t+\frac9{32}\right)e^{-4t}.  \tag{3}
\]

Thus \(V_{n,i}^2\) is bounded by an outward enclosure of
\(F_n(a_i)-F_n(b_i)\). The \(V_{2,i}\) caps are taken from the
frozen exact \(w_2\) profile data in the passing unit-W certificate.
No estimate of \(w_3\) is substituted for an earlier iterate.

## 3. Three finite recurrences prove entry

Take the certified rational upper constant
\(\bar c=707106781187/2000000000000>c\), and define nonnegative
vectors recursively by

\[
 E_{0,i}=s_i,
\]
\[
 E_{n+1,i}\ge s_i+
 \bar c\sum_{j\le i}\kappa_{ij}
 \left(2V_{n,j}E_{n,j}+E_{n,j}^2\right),
 \qquad n=0,1,2.                                  \tag{4}
\]

All rounding in (4) is outward. Subtracting the two exact Picard
formulas gives

\[
 v_h^{[n+1]}-w_{n+1}
 =e^{-tA}h-2\mathcal Q(w_n,v_h^{[n]}-w_n)
  -\mathcal Q(v_h^{[n]}-w_n,v_h^{[n]}-w_n).
\]

The causal block bilinear estimate and induction therefore prove

\[
 \|v_h^{[n]}-w_n\|_{\mathcal X(I_i)}\le E_{n,i}
 \quad(0\le n\le3).                               \tag{5}
\]

The independently certified perturbed map \(T_h\) preserves the
product ball centered at the unperturbed \(w_3\), with block radii
\(\rho_i\), for every \(A_{-1}(h)\le1/10000\). Its certificate
retains the full initial heat correction in its defect. The fresh
finite-prefix calculation verifies the exact strict inequalities

\[
 \boxed{E_{3,i}<\rho_i\quad\text{for all 17 blocks}.}       \tag{6}
\]

For orientation, the largest ratio is
\(\max_i E_{3,i}/\rho_i=0.9127337584808773\ldots<1\).
The comparison itself is between exact rational numbers. Equations
(5)--(6) place the intended heat-seeded third iterate in the invariant
ball; they do not merely show existence of some starting point there.

## 4. Ordinary convergence and an honest error bound

The perturbed certificate supplies a nonnegative lower triangular
Lipschitz matrix \(M\) with

\[
 M_{ij}=2\bar c\kappa_{ij}(V_{3,j}+\rho_j)\ (j<i),
 \qquad M_{ii}=2\bar c(V_{3,i}+\rho_i),
\]

and \(\max_iM_{ii}<3/5\). Positive rational weights for any
\(q>\max_iM_{ii}\) give an equivalent finite-block norm in which
\(T_h\) is a strict contraction. In particular take \(q=3/4\).
The invariant product ball is closed in the complete space
\(\mathcal X\). The mild map supplies the continuous critical trace
with initial datum \(W+h\); convergence of its bilinear expressions
also preserves that trace. Thus Banach's theorem, together with (6),
proves convergence of the actual heat-seeded sequence to the fixed
point \(u_h\).

An error statement which keeps the finite-prefix contribution is

\[
 \boxed{
 \|v_h^{[n]}-u_h\|_{\mathcal X(I_i)}
 \le\left[M^{n-3}(E_3+\rho)\right]_i,
 \qquad n\ge3.}                                   \tag{7}
\]

Consequently the global \(\mathcal X\)-error is bounded by the
Euclidean norm of the vector on the right. The initial fixed-point
error in (7) is \(E_3+\rho\), rather than just \(\rho\); using the
unperturbed-centered error bound without this extra term would omit
the perturbed prefix. A rational geometric prefactor can be obtained
from any positive vector \(p\ge E_3+\rho\) satisfying
\(Mp\le(3/4)p\).

The companion fresh certificate
`perturbed_stability_rate_static_result.json` constructs such positive
majorants for several choices of rate. In particular, it proves the
simple uniform headline bound

\[
 \boxed{\|v_h^{[n]}-u_h\|_{\mathcal X}
       \le6(9/10)^{n-3},\qquad n\ge3.}             \tag{8}
\]

Two additional valid bounds are \(111(4/5)^{n-3}\) and
\(961(3/4)^{n-3}\); their smaller geometric factors have larger
prefactors. Each is proved from \(p\ge E_3+\rho\) and
\(Mp\le q p\), with a rational enclosure of \(\|p\|_{\ell^2}\).
They improve the valid but coarse all-ones-resolvent prefactor
recorded in the prefix JSON, without changing the finite-prefix data.

The same companion certificate gives
\(\|u_h-u_g\|_{\mathcal X}\le4757A_{-1}(h-g)\)
for two perturbations satisfying the stated radius. This follows from
the common block ball and the positive triangular resolvent
\((I-M)^{-1}\) applied to the initial heat-source caps; it does not
require the perturbed solutions to share the same initial datum.

The earlier result for \(|a|\le1/2\), \(A_{-1}(h)\le1/100\),
used a different small-reference contraction. The present claim is an
explicit neighborhood of the unit-amplitude datum and concerns the
true full-data heat seed. Neither statement establishes ordinary
Picard convergence on the full prior regularity tube
\(A_{-1}(h)\le1/20\).

## Verification record

`perturbed_picard_prefix_static.py` is new standalone code using only
Python's standard library. Its 493 exact assertions check the input
certificate linkage, analytic endpoint enclosures, all three finite
recurrences, all 17 entry comparisons, and a positive weighted-matrix
contraction construction. It writes the full \(E_0,E_1,E_2,E_3\)
vectors and \(E_3+\rho\) to
`perturbed_picard_prefix_static_result.json`, with status
`PASS_PERTURBED_TRUE_HEAT_SEED_FINITE_PREFIX_ENTRY`.

The frozen source `perturbation_block_gate_result.json` has SHA-256
`c789c350a041c8e9eba01ccea50f1839fbe2a678439f924f6c13e4f29f2ca461`.
The prefix result is frozen at SHA-256
`0b6be7a25f469189036951ec5bfcea2a1cf62c35d8a2bbe50ee838ea2645e77e`.
The companion rate and stability checker records 441 passing exact
assertions and is bound to those same two source hashes.
Earlier files are unchanged. The finite-prefix check is a static
analytic certificate, not numerical time integration.
