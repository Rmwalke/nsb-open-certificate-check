# Independent audit: perturbing the unit-W Picard certificate

This note concerns ordinary global Picard iteration for initial data
\(u_0=W+h\), where \(W=(\sin z,\sin x,\sin y)\) and \(h\) is real,
mean-zero and divergence-free. To compare directly with the previously
established global regularity class, take \(h\) smooth and measure it by
\(\chi(h)=\sum_{k\ne0}|k|^{-1}|\widehat h(k)|\).
The existing unit-W block certificate is an input, not recomputed here.

The independent audit passes for \(\chi(h)\le1/10000\). The actual
heat-seeded sequence \(v_0=S(W+h)\), \(v_{n+1}=S(W+h)-Q(v_n,v_n)\)
converges and satisfies
\(\|v_n-u_h\|_{\mathcal X}\le6(9/10)^{n-3}\) for every \(n\ge3\).
This neighborhood lies inside the already established regularity class
\(\chi(h)\le1/20\). The result is a new iteration-convergence certificate,
not an enlargement of that regularity class.

## 1. Full heat contribution of the initial perturbation

Write \(\mathcal X(I)=L^2(I;A_0)\), \(Q=HB_s\),
\(c=1/(2\sqrt2)\), and \(T_h(v)=S(W+h)-Q(v,v)\).
For an output block \(I=[a,a+\ell]\), Fourier-mode Minkowski gives

\[
\|Sh\|_{\mathcal X(I)}
\le\frac{\chi(h)}{\sqrt2}\,
\sup_{\lambda\ge1}e^{-\lambda a}
                \sqrt{1-e^{-2\lambda\ell}}.
\]

The last block has \(\ell=\infty\). This is the same dimensionless kernel
factor as in the frozen block certificate, with gap \(a\). Its coefficient
is \(1/\sqrt2\), twice the nonlinear coefficient \(c\). Use a rigorous
upper bound \(H_i\) for \(\|Sh\|_{\mathcal X(I_i)}\), retaining the entire
initial perturbation at every output time.

For the first block, the supremum is one for every positive length: a
perturbation may be concentrated at arbitrarily high frequencies. A
uniform factor tending to zero with the first-block length would require
an additional frequency or stronger-norm restriction. For any fixed
\(h\) of finite \(\chi\), its actual first-block norm does tend to zero as
the block length tends to zero; that pointwise assertion does not give a
uniform small factor over a critical-norm ball.

## 2. An invariant rectangle around the unperturbed third reference

Keep the same exact unperturbed \(w_3\), with complete heat defect \(z_4\).
Let \(V_i\ge\|w_3\|_{\mathcal X(I_i)}\),
\(D_i\ge\|z_4\|_{\mathcal X(I_i)}\), and let \(\kappa_{ij}\) be the
frozen causal block-kernel upper bounds, with \(\kappa_{ii}=1\).
Then

\[
T_h(w_3+e)-w_3=Sh-z_4-2Q(w_3,e)-Q(e,e).
\]

Given a rational \(\bar c\ge c\), positive radii \(\rho_i\) give an
invariant rectangle if

\[
D_i+H_i+
\bar c\sum_{j<i}\kappa_{ij}(2V_j\rho_j+\rho_j^2)
+2\bar cV_i\rho_i+\bar c\rho_i^2\le\rho_i
\quad\text{for every }i.
\]

For two points in that rectangle, their blockwise image differences are
bounded by the lower triangular matrix

\[
A_{ii}=2\bar c(V_i+\rho_i),\qquad
A_{ij}=2\bar c\kappa_{ij}(V_j+\rho_j)\ (j<i).
\]

If every diagonal is strictly less than one, positive weights make this a
contraction. The product of closed block balls is complete; the map is the
same full-history global map on every block. The heat perturbation is an
additive forcing and cancels in the difference of two map values, so it
does not directly appear in \(A\). Its effect on the permitted radii must
still be included.

This argument by itself only proves that the rectangle contains a unique
fixed point and that seeds already inside it converge. It does not yet
place the requested perturbed heat seed, or one of its iterates, inside.

## 3. Certifying the actual perturbed heat-seeded finite prefix

Let \(w_0=SW\), \(w_{n+1}=T_0(w_n)\), and let
\(v_0=S(W+h)\), \(v_{n+1}=T_h(v_n)\). The previously constructed exact
references are precisely \(w_1,w_2,w_3\). Put \(\delta_n=v_n-w_n\).
Then

\[
\delta_0=Sh,\qquad
\delta_{n+1}=Sh-2Q(w_n,\delta_n)-Q(\delta_n,\delta_n).
\]

For upper bounds \(V_{n,j}\ge\|w_n\|_{\mathcal X(I_j)}\), the recursion

\[
E_{0,i}=H_i,\qquad
E_{n+1,i}=H_i+
\bar c\sum_{j\le i}\kappa_{ij}
                  (2V_{n,j}E_{n,j}+E_{n,j}^2)
\]

is valid for \(n=0,1,2\). One may replace this formula by any independently
justified sharper upper bound. If \(E_{3,i}\le\rho_i\) for every block,
then the actual third perturbed iterate belongs to the invariant rectangle
and all its later iterates converge to its fixed point. There are no
restarts and no discarded initial tail.

If \(E_3\le\rho\), the error of this seed relative to the fixed point is
bounded initially by \(\rho+E_3\), not automatically by \(\rho\).
Consequently a valid matrix error estimate is

\[
\|v_n-u\|_{\mathcal X}
\le\|A^{n-3}(\rho+E_3)\|_{\ell^2},\qquad n\ge3.
\]

The simpler vector \(2\rho\) is also safe. Any old scalar rate constant
must therefore be rechecked for the new radii and this new initial error
vector. A scalar rate \(q\) must exceed all of the new diagonal entries;
the unit-W weights or constant are not automatically valid unchanged.

## 4. Trace, identification, and scope

The fixed-point equation is \(u=S(W+h)-Q(u,u)\), so it has critical initial
trace \(W+h\). The nonlinear heat term tends to zero in the critical norm
at the initial endpoint, since its bound is controlled by the local
\(\mathcal X\)-norm squared. Global \(\mathcal X\) uniqueness follows by
successively choosing intervals where two solutions have sufficiently
small local norms; equality of their pasts removes the earlier history
from their difference equation. For smooth perturbations in the already
regular class, this identifies the limit with that regular solution.

The certificate is uniform over the specified critical-norm ball; it is
not a claim about unrestricted perturbations, arbitrary Picard seeds,
arbitrary Navier-Stokes initial data, or convergence following merely from
global regularity. Finite Fourier support is not assumed for \(h\).

## 5. Concrete data audit

The new data-only program `perturbed_picard_adversarial_static.py` returned
`PASS_PERTURBED_TRUE_HEAT_SEED_INDEPENDENT_AUDIT`. It imports and executes no
earlier research program. It checks the frozen input hashes, independently
encloses exponentials by positive Taylor sums with rational remainders,
and checks the following exact rational statements:

- All 17 initial heat caps are valid, including the first factor one.
- All 17 enlarged invariant-rectangle inequalities have strictly positive
  slack. The smallest slack is approximately \(1.15\times10^{-19}\);
  exact rational arithmetic establishes its positive sign.
- Every new matrix entry has the required causal formula and every
  diagonal is below \(3/5\). The largest diagonal is approximately
  \(0.5781572417201134\).
- All 51 reference norm caps for \(w_0,w_1,w_2\) and all 51 finite-prefix
  update caps hold. The actual third perturbed iterate lies strictly
  inside every block ball; the maximum ratio \(E_{3,i}/\rho_i\) is
  approximately \(0.9127337584808773\).
- The exact solution-budget sum is less than nine (approximately
  \(8.340335743815434\)).

For the first two reference caps, the independently checked squared-norm
tails are

\[
\int_t^\infty A_0(w_0)^2=\frac92e^{-2t},
\]
\[
\int_t^\infty A_0(w_1)^2=
\frac92e^{-2t}+(6t+2)e^{-3t}
+\left(\frac94t^2+\frac98t+\frac9{32}\right)e^{-4t}.
\]

The \(w_2\) caps are tied to the already audited frozen exact profile.
Interval moments are obtained by subtracting the two endpoint tails with
the correct outward rounding. No time quadrature or evolution was used.

The principal frozen input hashes are:

- `perturbation_block_gate_result.json`:
  `c789c350a041c8e9eba01ccea50f1839fbe2a678439f924f6c13e4f29f2ca461`.
- `perturbed_picard_prefix_static_result.json`:
  `0b6be7a25f469189036951ec5bfcea2a1cf62c35d8a2bbe50ee838ea2645e77e`.
- `perturbed_stability_rate_static_result.json`:
  `76536c8207a28198339dc502fe1f023869a00e140fb2d09f4d1b13e6a3549499`.

The complete bindings and check results are saved in
`perturbed_picard_adversarial_static_result.json`.

## 6. Quantitative convergence and solution-map stability

Set \(E=\rho+E_3\). For any \(q\) above the diagonal of \(A\), the recursion

\[
p_i=\max\left\{E_i,
 \frac{\sum_{j<i}A_{ij}p_j}{q-A_{ii}}\right\}
\]

produces a positive vector with \(p\ge E\) and \(Ap\le qp\). Thus
\(A^jE\le q^jp\), and the global \(\mathcal X\) prefactor is at most
\(\|p\|_{\ell^2}\). The independent checker verifies every component of
these inequalities and the exact squared prefactor bounds. For \(n\ge3\),
the following three simultaneous bounds are certified:

| Rate | Valid global error bound |
|---|---|
| \(3/4\) | \(961(3/4)^{n-3}\) |
| \(4/5\) | \(111(4/5)^{n-3}\) |
| \(9/10\) | \(6(9/10)^{n-3}\) |

These bounds use the complete initial error \(\rho+E_3\). The old unit-W
scalar constant is not reused without a new check.

For two perturbations \(h,g\) in the certified ball, both fixed points lie
in the same rectangle. Let \(z_i=\|u_h-u_g\|_{\mathcal X(I_i)}\), and let
\(s_i\) denote the heat cap per unit \(\chi\). Their difference satisfies
\(z\le s\chi(h-g)+Az\). As \(A\) is lower triangular with diagonal below
one, successive substitution gives

\[
z_i\le b_i\chi(h-g),\qquad
b_i=\frac{s_i+\sum_{j<i}A_{ij}b_j}{1-A_{ii}}.
\]

The data-only checker verifies \((I-A)b=s\), positivity, and
\(\sum_i b_i^2\le4757^2\). Consequently

\[
\|u_h-u_g\|_{\mathcal X}\le4757\,\chi(h-g).
\]

This bound requires each of \(h\) and \(g\) to lie in the certified ball.
It does not require their difference to lie in that same ball: the heat
estimate is homogeneous and the difference inequality is linear in its
critical initial norm.

## 7. Budget and the uniform Galerkin NS-B estimate

The invariant rectangle gives
\(\int_0^\infty A_0(u_h)^2\le\sum_i(V_i+\rho_i)^2<9\).
For every Galerkin cutoff \(K\ge\sqrt{24}\), the same proof applies to the
projected initial data: \(w_3\) is supported in squared frequencies at
most 22, projections contract the Fourier norms, and
\(\chi(P_Kh)\le\chi(h)\). The projected residual still contains its full
history. Thus the Fourier budget below nine holds uniformly over these
large cutoffs.

Put \(E_0=\|W+h\|_2^2\) and
\(X_0=\|\nabla(W+h)\|_2^2\). The previously proved energy/enstrophy
estimate yields

\[
\int_0^\infty\|\nabla u_K(t)\|_2^4\,dt
\le\frac{E_0X_0}{2}e^9
\quad (K\ge\sqrt{24}).
\]

For \(K<\sqrt{24}\), energy alone gives
\(\int\|\nabla u_K\|_2^4\le K^2E_0^2/2<12E_0^2\).
Mean-zero Poincare gives \(E_0\le X_0\), while
\(e^9>1+9+81/2>24\), so the same final bound absorbs these small cutoffs.
The NS-B estimate is therefore uniform over all cutoffs. Its prefactor
depends on the smooth initial datum through \(E_0X_0\); a critical-norm
ball does not impose a uniform bound on that higher-norm product.

No substantive flaw was found in the stated \(\chi\le1/10000\) scope.
Failure of the tested larger-radius majorants is not evidence of Picard
divergence. Extending convergence to the older full regularity ball
\(\chi\le1/20\) remains uncertified by this calculation.
