# Explicit rates and solution-map stability in the certified perturbation tube

Work in the normalized periodic setting of the companion proof. Let
\[
\chi(h)=A_{-1}(h),\qquad \varepsilon=10^{-4},\qquad
\mathcal X=L^2(0,\infty;A_0).
\]
All perturbations are smooth, real, mean-zero and divergence-free. The frozen block certificate constructs a single closed rectangle
\[
\mathcal B=\{f:\|f-w_3\|_{\mathcal X(I_i)}\le R_i,\ 1\le i\le17\}
\]
invariant under every \(T_h(f)=S(W+h)-Q(f,f)\) with \(\chi(h)\le\varepsilon\).
Its nonnegative lower-triangular Lipschitz matrix is
\[
A_{ij}=2\bar c\,k_{ij}(V_j+R_j)\quad(j\le i),\qquad A_{ij}=0\quad(j>i).
\]
The independent finite-prefix calculation bounds the third heat-seeded iterate:
\[
\|v_h^{(3)}-w_3\|_{\mathcal X(I_i)}\le E_{3,i}<R_i,
\quad v_h^{(0)}=S(W+h),\quad v_h^{(n+1)}=T_h(v_h^{(n)}).
\]
Its initial error relative to the fixed point is therefore bounded by
\(E_i=E_{3,i}+R_i\). The \(R_i\) term cannot be omitted.

## 1. Smaller rigorous scalar prefactors

For a fixed \(q\in(\max_i A_{ii},1)\), define recursively
\[
p_i=\max\left\{E_i,\,
       \frac{\sum_{j<i}A_{ij}p_j}{q-A_{ii}}\right\}.
\]
Then \(p_i\ge E_i>0\) and \(Ap\le qp\), directly by construction. This is the componentwise least vector satisfying these two inequalities: any other such vector dominates it by induction in \(i\). This observation concerns this particular certificate system, not optimal convergence rates of the actual equation.

If \(e_j\) is the vector of block errors after \(j\) further iterations, positivity gives
\[
e_j\le A^jE\le A^jp\le q^jp.
\]
Because the time blocks partition the half-line,
\[
\|v_h^{(n)}-u_h\|_{\mathcal X}
\le\|p\|_{\ell^2}q^{n-3}\qquad(n\ge3).
\]
All operands and comparisons in the new checker are rational. The following integer prefactors are independently valid:

| \(q\) | Certified integer prefactor | Exact squared norm is below |
|---|---:|---:|
| \(3/4\) | 961 | \(961^2\) |
| \(4/5\) | 111 | \(111^2\) |
| \(9/10\) | 6 | \(6^2\) |

Thus a simple uniform statement is
\[
\boxed{\|v_h^{(n)}-u_h\|_{\mathcal X}
       \le6(9/10)^{n-3}\quad(n\ge3).}
\]
One may use the minimum of the three displayed geometric bounds. The exact componentwise estimate \(A^{n-3}(E_3+R)\) is also available. The factors are upper bounds from the finite matrix; they do not assert an observed numerical convergence speed. No \(H^1\) convergence rate is inferred.

## 2. Uniform dependence on the initial datum

Let \(\chi(h),\chi(g)\le\varepsilon\), and write \(u_h,u_g\) for the two fixed points in the same rectangle. The initial heat factors \(s_i\) satisfy
\[
\|S(h-g)\|_{\mathcal X(I_i)}\le s_i\chi(h-g).
\]
The difference of the fixed-point equations consequently gives, componentwise,
\[
z\le s\,\chi(h-g)+Az,\qquad
z_i=\|u_h-u_g\|_{\mathcal X(I_i)}.
\]
Define the positive vector
\[
b_i=\frac{s_i+\sum_{j<i}A_{ij}b_j}{1-A_{ii}}.
\]
Forward substitution, using \(A_{ii}<1\), proves \(z_i\le b_i\chi(h-g)\).
Equivalently, \(b=(I-A)^{-1}s\); its entries are nonnegative without any appeal to an unsigned inverse estimate.

The exact calculation gives
\[
\sum_i b_i^2<4757^2,
\]
and hence
\[
\boxed{\|u_h-u_g\|_{\mathcal X}\le4757\,A_{-1}(h-g).}
\]
Both endpoints must belong to the certified tube. The perturbation difference itself may have norm up to \(2\varepsilon\). The constant is conservative; it is not a lower bound on amplification or evidence of instability.

## 3. Uniform NS-B bound and verification

The common rectangle has
\[
\sum_i(V_i+R_i)^2=8.340335743815434\ldots<9.
\]
For every Fourier Galerkin cutoff containing \(w_3\), in particular \(K\ge\sqrt{24}\), projection preserves the certificate, so
\(\int A_0(u_{h,K})^2<9\).
Energy and enstrophy then give
\[
\int_0^\infty\|\nabla u_{h,K}\|_2^4dt\le\frac{E_0X_0}{2}e^9.
\]
For \(K<\sqrt{24}\), energy and mean-zero Poincaré give instead
\[
\int_0^\infty\|\nabla u_{h,K}\|_2^4dt
\le\frac{K^2E_0^2}{2}\le12E_0X_0
<\frac{E_0X_0}{2}e^9.
\]
Here \(E_0=\|W+h\|_2^2\) and \(X_0=\|\nabla(W+h)\|_2^2\) are the full initial-data quantities, which need not be close to their values at \(W\) under the critical-norm assumption. The budget below nine is not claimed for every small cutoff.

The fresh standard-library checker is **perturbed_stability_rate_static.py**, with exact rational output in **perturbed_stability_rate_static_result.json**. It passed 441 finite assertions, including source binding, actual prefix entry, all matrix formulas, each majorant inequality, both forward-substitution identities and every stated integer constant. It reads the new gate and prefix JSON files as data. No historical program is executed or imported; no evolution equation or quadrature is evaluated.

The older regularity tube \(\chi(h)\le1/20\) already contains this family. The new conclusions are ordinary heat-seeded iteration convergence, its quantitative rates, solution-map stability, and an improved analytic NS-B bound on the smaller tube. This is not an arbitrary-data theorem, external peer review, or formal proof-assistant verification.
