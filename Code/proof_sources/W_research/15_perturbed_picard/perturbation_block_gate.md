# A certified critical perturbation radius for residual-correction convergence

Private local continuation. Work on the normalized \(2\pi\)-periodic three-dimensional torus with viscosity one. Let
\[
W=(\sin z,\sin x,\sin y),\qquad
\chi(h)=A_{-1}(h)=\sum_{k\ne0}\frac{|\widehat h_k|}{|k|}.
\]
The perturbation \(h\) is real, smooth, mean-zero, and divergence-free. It need not have bounded Fourier support.

**Certified radius.** For every such \(h\) with
\[
\boxed{\chi(h)\le\varepsilon:=1/10000,}                         \tag{1}
\]
the ordinary mild Picard map for the complete datum \(W+h\),
\[
T_h(f)=e^{t\Delta}(W+h)-\mathcal Q(f,f),
\]
has a common contracting invariant block set centered at the frozen unit-data reference \(w_3\). Iteration from \(w_3\), or from \(w_3+e^{t\Delta}h\), converges in
\(\mathcal X=L^2(0,\infty;A_0)\) to the solution \(u_h\). In the second case every iterate has initial datum \(W+h\), and each update is literally the full residual heat correction.

The uniform solution bound is
\[
\boxed{\int_0^\infty A_0(u_h)^2dt<9.}                          \tag{2}
\]
For every solution Galerkin cutoff \(M\),
\[
\boxed{\int_0^\infty\|\nabla u_{h,M}(t)\|_2^4dt
 \le\frac{E_0X_0}{2}e^9,\quad
 E_0=\|W+h\|_2^2,\quad X_0=\|\nabla(W+h)\|_2^2.}               \tag{3}
\]
The separate finite-prefix certificate now also proves convergence of the standard heat-seeded sequence, starting from \(e^{t\Delta}(W+h)\). Section 3A explains its entry into this same invariant set; the initial heat seed itself is not assumed to lie in the set.

The older regularity tube already covered \(\chi(h)\le1/20\) at unit amplitude. Thus (1) is a new iteration-convergence certificate and a better analytic budget on a smaller family. It is not a new global-regularity class or a larger perturbation domain.

## 1. Frozen data and the complete initial contribution

Use the reference, grid, and full source-memory factors from the frozen certificate in *refined_W_block_gate_result.json*. Its center is
\[
v=w_3=w_2-\mathcal H r_2,\qquad v(0)=W,
\]
and its complete next heat defect is
\[
z_4=\mathcal H r_3,\qquad
r_3=-2B_s(w_2,z_3)+B_s(z_3,z_3),\quad z_3=\mathcal H r_2.
\]
Here
\[
\mathcal H r(t)=\int_0^t e^{(t-s)\Delta}r(s)ds,\qquad
\mathcal Q(a,b)=\mathcal H B_s(a,b),
\]
and \(B_s\) is the symmetric Leray-projected NS bilinear operator. All low residual modes and all earlier forcing times are retained.

The 17-block endpoints are
\[
0,\tfrac1{16},\tfrac18,\tfrac3{16},\tfrac14,\tfrac5{16},
\tfrac38,\tfrac7{16},\tfrac12,\tfrac58,\tfrac34,\tfrac78,
1,\tfrac54,\tfrac32,2,3,\infty.
\]
For \(I_i=[a_i,b_i]\), the frozen rational data supply
\[
V_i\ge\|v\|_{\mathcal X(I_i)},\qquad
d_i\ge\|z_4\|_{\mathcal X(I_i)}.
\]
The complete target-data center defect is
\[
T_h(v)-v=e^{t\Delta}h-z_4.                                    \tag{4}
\]
In particular, omitting \(e^{t\Delta}h\) or keeping the original datum \(W\) in the map would certify the wrong problem.

For \(\ell=b-a\), define
\[
\kappa(a,\ell)=\sup_{\lambda\ge1}
 e^{-\lambda a}\sqrt{1-e^{-2\lambda\ell}},
\]
with \(1-e^{-2\lambda\infty}=1\). Minkowski in the Fourier sum gives
\[
\begin{split}
\|e^{t\Delta}h\|_{\mathcal X(a,b)}
&\le\sum_{k\ne0}|\widehat h_k|
       \left(\int_a^b e^{-2|k|^2t}dt\right)^{1/2}\\
&=\frac1{\sqrt2}\sum_{k\ne0}
  \frac{|\widehat h_k|}{|k|}
 e^{-|k|^2a}\sqrt{1-e^{-2|k|^2\ell}}\\
&\le\frac{\chi(h)}{\sqrt2}\kappa(a,\ell).                      \tag{5}
\end{split}
\]
For the first block \(a=0\), the supremum is one, regardless of how short the block is. Arbitrarily high-frequency perturbations prevent a uniform square-root-of-duration improvement there.

The checker takes
\[
\bar c=\frac{707106781187}{2000000000000}>
c:=\frac1{2\sqrt2},
\qquad
s_i=2\bar c\,\overline{\kappa(a_i,b_i-a_i)}.                   \tag{6}
\]
Thus the full initial contribution is bounded by
\(\chi(h)s_i\), and (4) has block norm at most \(d_i+\varepsilon s_i\). The maximum \(s_1=2\bar c\) is an upper enclosure of \(1/\sqrt2\); the later factors are strictly smaller.

For \(a,\ell>0\), the unconstrained maximum over \(\lambda\ge0\) of the square of the factor is
\[
\frac{\ell}{a+\ell}
 \left(\frac a{a+\ell}\right)^{a/\ell}.
\]
If the maximizing \(\lambda_*=\log(1+\ell/a)/(2\ell)\) is at most one, the constrained maximum is attained at one. The new checker encloses these expressions by rational Taylor, logarithm-series, and square-root bounds. The infinite tail uses \(\kappa(a,\infty)=e^{-a}\).

## 2. Common invariant radii for the entire perturbation family

The established causal bilinear bound is
\[
\|\mathcal Q(f,g)\|_{\mathcal X(I_i)}
 \le c\sum_{j\le i}k_{ij}
       \|f\|_{\mathcal X(I_j)}\|g\|_{\mathcal X(I_j)},          \tag{7}
\]
where the frozen rational \(k_{ij}\) include the full earlier source history. The diagonal is \(k_{ii}=1\). For \(j<i\), the factor bounds
\[
\sup_{\lambda\ge1}
 \sqrt{e^{-2\lambda(a_i-b_j)}
       (1-e^{-2\lambda(b_i-a_i)})}.
\]
Only the nonzero torus spectral gap \(\lambda\ge1\) is used. No lower frequency bound is assumed for the corrected residual.

For the single worst-case radius \(\varepsilon\) in (1), choose positive rational \(\rho_i\) satisfying
\[
\boxed{d_i+\varepsilon s_i+
 \bar c\sum_{j\le i}k_{ij}(2V_j\rho_j+\rho_j^2)\le\rho_i.}     \tag{8}
\]
These same radii then work for every \(\chi(h)\le\varepsilon\).

Indeed, the closed product set
\[
\mathcal B=\{f:\|f-v\|_{\mathcal X(I_i)}\le\rho_i
                  \text{ for every }i\}
\]
is invariant for \(T_h\), by (4), (5), and (7). On this set the vector of block distances is bounded by the nonnegative lower triangular matrix
\[
A_{ij}=2\bar c\,k_{ij}(V_j+\rho_j)\quad(j\le i),\qquad
A_{ij}=0\quad(j>i).                                          \tag{9}
\]
Every diagonal entry is below one, so an equivalent weighted block norm makes all these maps contractions. Explicitly, for any \(q\) between \(\max_iA_{ii}\) and one,
\[
p=(qI-A)^{-1}\mathbf1>0,\qquad Ap=qp-\mathbf1\le qp.
\]
The norm \(\max_i\|f\|_{\mathcal X(I_i)}/p_i\) is equivalent to \(\mathcal X\), since the partition has finitely many blocks. This is a contraction for each full global Picard map, not an interval-by-interval replacement of the evolution.

For iteration from \(v=w_3\),
\[
\big(\|T_h^n(v)-u_h\|_{\mathcal X(I_i)}\big)_i\le A^n\rho.
                                                                    \tag{10}
\]
The shifted center \(v+e^{t\Delta}h\) also lies in \(\mathcal B\), because (8) implies \(\varepsilon s_i\le\rho_i\). Its initial error relative to \(u_h\) is at most \(\rho_i+\varepsilon s_i\), so
\[
\big(\|T_h^n(v+e^{t\Delta}h)-u_h\|_{\mathcal X(I_i)}\big)_i
 \le A^n(\rho+\varepsilon s).                                \tag{11}
\]
This second implementation preserves the target initial datum from its first reference onward. If the unshifted center \(v\) is used, the map must explicitly be \(T_h\); ordinary residual subtraction from \(v\) alone would instead preserve \(W\).

## 3. The explicit passing radii

All values in the next table are decimal orientations only. The certified operands, including all source-memory factors, are exact rationals in *perturbation_block_gate_result.json*.

| Block | Interval | \(s_i\), the coefficient of \(\chi(h)\) | Enlarged \(\rho_i\) |
|---:|---|---:|---:|
| 1 | 0–1/16 | 0.707106781187 | 0.0001513601476 |
| 2 | 1/16–1/8 | 0.353553390594 | 0.0002615322810 |
| 3 | 1/8–3/16 | 0.272165526976 | 0.0005075400416 |
| 4 | 3/16–1/4 | 0.229639663386 | 0.0009973137978 |
| 5 | 1/4–5/16 | 0.202385770251 | 0.001862822387 |
| 6 | 5/16–3/8 | 0.183002179408 | 0.003239537389 |
| 7 | 3/8–7/16 | 0.168304455548 | 0.005247174213 |
| 8 | 7/16–1/2 | 0.156663633260 | 0.007971605592 |
| 9 | 1/2–5/8 | 0.201710986700 | 0.01837514495 |
| 10 | 5/8–3/4 | 0.178009320980 | 0.03400893713 |
| 11 | 3/4–7/8 | 0.157092674396 | 0.05347446314 |
| 12 | 7/8–1 | 0.138633798574 | 0.07516228487 |
| 13 | 1–5/4 | 0.163172124784 | 0.1394382766 |
| 14 | 5/4–3/2 | 0.127078578557 | 0.1994588479 |
| 15 | 3/2–2 | 0.125442077229 | 0.3090665577 |
| 16 | 2–3 | 0.088985635673 | 0.3816030484 |
| 17 | 3–infinity | 0.035204773658 | 0.2765245298 |

All 17 inequalities (8) pass. The minimum quadratic discriminant is greater than \(0.17795\), and
\[
\boxed{\max_i A_{ii}
 =0.5781572417201134\ldots<3/5.}                              \tag{12}
\]
The source-history inequalities have not discarded any prior time block.

The common invariant set satisfies
\[
\boxed{\sum_i(V_i+\rho_i)^2
 =8.340335743815434\ldots<9.}                                 \tag{13}
\]
This proves (2), uniformly over the entire family (1). The assertion does not bound \(\|h\|_{H^1}\) by its critical norm; high-frequency perturbations are permitted.

For an elementary explicit scalar rate, the checker sets \(q=3/4\) and uses the positive resolvent vector above. The resulting conservative bounds are
\[
\|T_h^n(w_3)-u_h\|_{\mathcal X}\le763(3/4)^n,
\]
\[
\|T_h^n(w_3+e^{t\Delta}h)-u_h\|_{\mathcal X}
 \le1119(3/4)^n.
\]
These constants are supplementary. The exact componentwise bounds (10) and (11) are usually more informative, and a separate finite-prefix/rate certificate may improve the scalar prefactors.

## 3A. Verified entry of the standard heat-seeded sequence

Let \(v_0^h=e^{t\Delta}(W+h)\), \(v_{n+1}^h=T_h(v_n^h)\), and let \(w_n\) be the unperturbed heat-seeded sequence. The independent fresh checker in *perturbed_picard_prefix_static.py* proves entry after three steps, using only exact endpoint moments and the same source-memory factors.

Specifically, let \(V_i^{(n)}\) bound \(\|w_n\|_{\mathcal X(I_i)}\) for \(n=0,1,2\). The first two reference amplitudes are exactly \(3e^{-t}\) and \(3e^{-t}+3te^{-2t}\), so their squared block norms have elementary exponential-polynomial antiderivatives. The third caps come from the frozen exact \(U_2\) certificate. Starting with \(E_i^{(0)}=\varepsilon s_i\), propagate
\[
E_i^{(n+1)}=\varepsilon s_i+
\bar c\sum_{j\le i}k_{ij}
 \left(2V_j^{(n)}E_j^{(n)}+(E_j^{(n)})^2\right),
\qquad n=0,1,2,
\]
using outward rational upper bounds. Formula (7) proves
\(\|v_n^h-w_n\|_{\mathcal X(I_i)}\le E_i^{(n)}\). All earlier forcing is present in each sum.

The certified result is
\[
\max_i\frac{E_i^{(3)}}{\rho_i}
=0.91273375848\ldots<1.
\]
Thus \(v_3^h\in\mathcal B\), and the ordinary global heat-seeded iteration converges. Its precise componentwise error bound is
\[
\boxed{\big(\|v_n^h-u_h\|_{\mathcal X(I_i)}\big)_i
 \le A^{n-3}(E^{(3)}+\rho),\qquad n\ge3.}
\]
The addition of \(E^{(3)}\) is necessary: the entry reference is the perturbed third iterate, rather than the unperturbed center \(w_3\). The exact prefix vectors and their positive margins are stored in *perturbed_picard_prefix_static_result.json*. This finite algebraic entry argument supplies the ordinary heat-seed conclusion without changing the iteration or integrating an evolution equation numerically.

## 4. Uniform Galerkin and NS-B consequence

The reference \(v=w_3\) has support \(|k|^2\le24\). For every cutoff \(M\ge\sqrt{24}\), \(P_Mv=v\), and
\[
T_{h,M}(f)=e^{t\Delta}P_M(W+h)-P_M\mathcal Q(f,f).
\]
Its center defect is \(P_M(e^{t\Delta}h-z_4)\). Projection contracts every Fourier-amplitude bound, including \(\chi(P_Mh)\le\chi(h)\). The same invariant set therefore constructs the Galerkin fixed point and gives
\[
\int_0^\infty A_0(u_{h,M})^2dt<9
\qquad(M\ge\sqrt{24}).                                      \tag{14}
\]
The construction is noncircular: it does not assume the full solution before obtaining the uniform estimates.

For \(X_M=\|\nabla u_{h,M}\|_2^2\) and
\(Y_M=\|\Delta u_{h,M}\|_2^2\), the established \(H^1\) inequality is
\[
X_M'+Y_M\le A_0(u_{h,M})^2X_M.
\]
Thus \(\sup_tX_M\le X_0e^9\), while the energy identity yields
\(\int_0^\infty X_M\le E_0/2\). Multiplying proves (3) for the cutoffs covered by (14).

For \(M<\sqrt{24}\), no Fourier-amplitude budget from the reference is asserted. Instead,
\[
\sup_tX_M\le24E_0,\qquad
\int_0^\infty X_M\le E_0/2,
\]
so
\[
\int_0^\infty X_M^2\le12E_0^2\le12E_0X_0
 <\frac{E_0X_0}{2}e^9.
\]
The middle step uses the mean-zero Poincaré inequality \(E_0\le X_0\). No smallness bound on the full \(E_0\) or \(X_0\) is required from the critical perturbation norm. This completes the all-cutoff claim.

The uniform strong estimates give global smooth continuation for each smooth datum in the family and identify that solution with the contracting mild fixed point.

## 5. Bounded search and what failures mean

A finite declared trial set tested the same grid, reference, source-memory factors, and initial heat profile. It did not search over new grids or optimize indefinitely.

| Trial \(\varepsilon\) | Declared majorant result | First failed block |
|---:|---|---:|
| 0 | Pass | — |
| 1/100000 | Pass | — |
| 1/50000 | Pass | — |
| 1/10000 | Pass | — |
| 1/5000 | Fail | 16 |
| 1/2000 | Fail | 13 |
| 1/1000 | Fail | 12 |
| 1/500 | Fail | 11 |
| 1/200 | Fail | 9 |
| 1/100 | Fail | 8 |
| 1/20 | Fail | 4 |

The selected value \(1/10000\) has a substantial diagonal and discriminant margin. It is not claimed to be the largest value the grid could certify. For each failed trial, lower bounds on the preceding minimal admissible radii lead to a negative upper discriminant in the indicated row. This proves failure of the specified positive majorant system, not failure of the actual iteration, loss of regularity, or impossibility of another certificate.

In particular, the failure at \(1/20\) does not contradict the older regularity result. Existence and convergence of one chosen iteration are different claims, and these componentwise upper bounds are conservative.

## 6. Reproducibility and frozen source binding

The new standalone standard-library program *perturbation_block_gate.py* reads only the frozen prior certificate as data. It does not execute or import historical research scripts. It first rechecks the imported baseline invariant-set algebra, then constructs the new initial heat caps, propagates exact rational lower and upper radius enclosures, and checks each final invariant-set and contraction inequality. All exponential and logarithmic enclosures are analytic rational bounds, and all square roots use integer-square-root enclosures. Decimal values have no role in acceptance.

| Artifact | SHA-256 |
|---|---|
| refined_W_block_gate_result.json, input | 17b578bd37cb971cfbd21a58749bc95a94d62e126fcced0fe21bb3b85a3420b0 |
| perturbation_block_gate_result.json, output | c789c350a041c8e9eba01ccea50f1839fbe2a678439f924f6c13e4f29f2ca461 |

The output status is **PASS_PERTURBED_W3_BLOCK_CONVERGENCE**. It exposes the selected radius, reference and defect caps, initial heat coefficients, full source-memory factors, passing rational radii, contraction matrix, and solution budget for independent reuse.

No numerical PDE/ODE evolution, quadrature, old-script execution/import, external search, upload, or change to older-stage artifacts was used. This is a proof supported by finite exact arithmetic, not verification inside a formal proof assistant.
