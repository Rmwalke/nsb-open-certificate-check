# Causal high-frequency residual memory with retained covariance

Private written NS-B continuation audit. No earlier research program, import, numerical integration, external source, upload, or original-file modification is used. This file proves estimates for the complete residual of an exact spherical Galerkin reference. It does not assume that the full PDE error remains at high frequencies.

Use viscosity one on the normalized 2π torus. All fields are real, mean-zero, divergence-free, and Fourier coefficients use the Euclidean complex-vector norm. Write Λ=(−Δ)^(1/2), A_j(w)=Σ_{k≠0}|k|^j|ŵ_k|, and

\[
\mathcal X(I)=L^2(I;A_0),\qquad
\mathfrak h(w)=\|e^{t\Delta}w\|_{\mathcal X(0,\infty)},\qquad
c=\frac1{2\sqrt2}.
\]

## 1. Full residual support and its retained input band

Let v=v_N solve the exact Galerkin equation with spherical support |k|≤N, N≥1. Set

\[
r=(I-P_N)B(v,v),\quad R=A_{-1}(r),\quad
U=A_0(v),\quad Y=\|\Lambda^2v\|_2^2,\quad Z=\|\Lambda^3v\|_2^2.
\]

Every residual output lies in N<|k|≤2N. Its inputs still lie inside the retained cutoff. An exterior output requires at least one retained input in N/2<|p|≤N, so the symmetric Leray estimate gives

\[
R\le U\,U_H-\tfrac12U_H^2\le U\,U_H
\le\sqrt{S_b}\,U\sqrt Z,\qquad
S_b=\sum_{N/2<|p|\le N}|p|^{-6}.                 \tag{1}
\]

Thus the band is retained, rather than replaced by the empty set of reference inputs above N. All statements below include the full exterior output shell. The safe high-output spectral gap λ=N² is used; the exact lattice minimum of |k|² among possible residual outputs can replace it by a larger value.

## 2. Exact causal damping for the defect alone

For 0≤τ<T define the recent forcing defect

\[
z_{\tau}(t)=\int_\tau^t e^{(t-s)\Delta}r(s)\,ds,
\quad \tau\le t\le T,
\]

and its scalar causal budget

\[
C_{\tau,T}=\int_\tau^T e^{-\lambda(T-s)}R(s)\,ds.
\]

Modewise heat decay, positivity after taking absolute values, and Minkowski give

\[
\boxed{A_{-1}(z_\tau(T))\le C_{\tau,T},\qquad
       \mathfrak h(z_\tau(T))\le C_{\tau,T}/\sqrt2.}             \tag{2}
\]

The second inequality follows from the exact future-time kernel norm

\[
\bigl\|e^{-|k|^2(T-s+t)}\bigr\|_{L^2(t>0)}
=\frac{e^{-|k|^2(T-s)}}{\sqrt2|k|}.
\]

For the total residual defect z₀, old forcing is not erased. It obeys

\[
A_{-1}(z_0(T))
\le e^{-\lambda(T-\tau)}A_{-1}(z_0(\tau))+C_{\tau,T},           \tag{3}
\]

with the analogous heat-norm bound

\[
\mathfrak h(z_0(T))
\le e^{-\lambda(T-\tau)}\mathfrak h(z_0(\tau))
   +C_{\tau,T}/\sqrt2.                                        \tag{4}
\]

If a forcing chunk ends at T₀, its subsequent contribution evolves by the heat equation. Both norms in (2) then decay by exp(−λ(t−T₀)). Its remaining solution-space norm after T₀+L is bounded by exp(−λL) times its endpoint heat norm. This is a precise short-memory theorem for the forcing defect.

## 3. A quadratic memory budget gains one power of N

Extend the forcing by zero outside I=[τ,T], and include the entire future heat response of that truncated forcing. Put f_k(s)=|r̂_k(s)|/|k| and λ_k=|k|², so Σ_k f_k=R and λ_k>λ=N². The positive pointwise majorant of its A₀ norm is

\[
\sum_k |k|\int_\tau^{\min(t,T)}e^{-\lambda_k(t-s)}f_k(s)\,ds.
\]

Square this nonnegative expression and use Tonelli. For fixed forcing times s≤r, the exact overlap integral of the two causal kernels is

\[
\int_r^\infty |k||\ell|e^{-\lambda_k(t-s)-\lambda_\ell(t-r)}dt
={|k||\ell|\over\lambda_k+\lambda_\ell}
  e^{-\lambda_k(r-s)}
\le\tfrac12 e^{-\lambda(r-s)}.
\]

For r≤s the same expression uses λ_ℓ in the exponential and has the same bound. Summing the nonnegative f_k(s)f_ℓ(r) therefore proves the complete quadratic memory estimate

\[
\boxed{
\|z_\tau\|_{\mathcal X(\tau,T)}^2
\le\|z_\tau^{\rm full\ future}\|_{\mathcal X(\tau,\infty)}^2
\le\mathcal M_{\tau,T}:=
\frac12\int_\tau^T\!\int_\tau^T
 e^{-\lambda|s-r|}R(s)R(r)\,ds\,dr.}                         \tag{5}
\]

Only the lower output gap is needed in this proof; the upper support bound 2N is unnecessary. The factor 1/2 comes from 2|k||ℓ|≤|k|²+|ℓ|². No reversed Minkowski inequality or singular kernel supremum is used.

The two-sided kernel exp(−λ|t|) has L¹ norm 2/λ. Scalar Young and Cauchy–Schwarz thus give

\[
\boxed{\|z_\tau\|_{\mathcal X(\tau,T)}
\le\sqrt{\mathcal M_{\tau,T}}
\le {1\over N}\|R\|_{L^2(\tau,T)}.}                          \tag{6}
\]

Dropping the exponential instead gives the independent bound

\[
\boxed{\sqrt{\mathcal M_{\tau,T}}
\le {1\over\sqrt2}\|R\|_{L^1(\tau,T)}.}                     \tag{7}
\]

Thus the actual quadratic memory budget is at least as good as each of these two coarser bounds, and it retains the separation in time between forcing events. The quantities R and its norms are finite-reference quantities. Formula (6) does not assert an a priori cutoff-uniform L² bound for R.

There is a useful scalar representation. Let

\[
C(t)=\int_\tau^t e^{-\lambda(t-s)}R(s)ds,
\quad C'=R-\lambda C,\quad C(\tau)=0.
\]

Symmetry of the double integral and then integration of CC′ give

\[
\mathcal M_{\tau,T}
=\int_\tau^T R(t)C(t)dt
=\frac{C(T)^2}{2}+\lambda\int_\tau^T C(t)^2dt.
\]

These are exact identities for the scalar majorant. The term C(T)²/2 is part of the full-future bound; one cannot subtract the separate upper bound on the actual endpoint heat norm to infer a sharper prefix estimate without another argument.

## 4. Endpoint covariance with the entire retained H² energy

The newly proved incompressible pair inequality gives

\[
Y'+2Z\le2U\sqrt{YZ}.                                           \tag{8}
\]

Restart at an arbitrary τ. Put

\[
a_\tau(s)=\int_\tau^sU(t)^2dt,\quad
E_\tau(s)=e^{-a_\tau(s)/2}Y(s),\quad
W_\tau(s)=e^{-a_\tau(s)/2}Z(s),
\]

and

\[
q_\tau(s)=\sqrt{W_\tau(s)}-\tfrac12U(s)\sqrt{E_\tau(s)}.
\]

Completing the square in (8) gives the exact differential inequality

\[
\boxed{E_\tau'+2q_\tau^2\le0.}                              \tag{9}
\]

Consequently

\[
E_\tau(s)\le Y(\tau),\qquad
\int_\tau^Tq_\tau^2ds
\le\tfrac12\bigl[Y(\tau)-E_\tau(T)\bigr].                    \tag{10}
\]

The bracket is nonnegative. The scalar q_τ need not be nonnegative; the following use of Cauchy–Schwarz remains valid for either sign.

Define two forcing-time moments

\[
J_1(\tau,T)=\int_\tau^T e^{-\lambda(T-s)}U(s)^2\sqrt{Y(s)}\,ds,
\]

\[
J_2(\tau,T)=\int_\tau^T e^{-2\lambda(T-s)+a_\tau(s)/2}U(s)^2\,ds.
\]

The identity √Z=(U√Y)/2+exp(a_τ/4)q_τ and (1) now prove

\[
\boxed{C_{\tau,T}\le\sqrt{S_b}\left[
\frac{J_1(\tau,T)}2+
\sqrt{\frac{Y(\tau)-e^{-a_\tau(T)/2}Y(T)}2\,J_2(\tau,T)}
\right].}                                                     \tag{11}
\]

This bound retains the whole time profiles of U and Y, the final remaining weighted H² energy, and the high-frequency decay of every forcing time. It can be restarted on a recent window and combined with (3)–(4), without deleting earlier forcing. Large dissipation already spent cannot simultaneously be counted as large remaining weighted energy: that covariance is preserved by (10).

If scalar upper bounds M≥sup_[τ,T]U² and L≥sup_[τ,T](U²√Y) are available, then with Δ=T−τ and A=a_τ(T),

\[
J_1\le L\frac{1-e^{-\lambda\Delta}}\lambda,
\quad
J_2\le e^{A/2}M\frac{1-e^{-2\lambda\Delta}}{2\lambda}.          \tag{12}
\]

The exact moments in (11) are generally preferable. Formula (12) simply shows the explicit N⁻² memory scale when such local reference envelopes are known.

## 5. A useful stronger unweighted bound follows as a corollary

Set λ=0 in the preceding integral argument solely to bound ∫R; this is a deliberately weakened damping estimate. Let A=a_τ(T) and E_end=e^(−A/2)Y(T). Since E_τ≤Y(τ), direct integration in a_τ gives

\[
\int_\tau^T U^2\sqrt Y\,ds
\le4\sqrt{Y(\tau)}(e^{A/4}-1),
\quad
\int_\tau^T e^{a_\tau/2}U^2ds=2(e^{A/2}-1).
\]

Therefore

\[
\boxed{\int_\tau^T R(s)ds\le\sqrt{S_b}\left[
2\sqrt{Y(\tau)}(e^{A/4}-1)
+\sqrt{[Y(\tau)-E_{\rm end}](e^{A/2}-1)}\right].}              \tag{13}
\]

In particular,

\[
\int_\tau^T R\le3\sqrt{S_bY(\tau)}\,e^{A/4}.                 \tag{14}
\]

At τ=0 and T→∞, (13)–(14) and (7) improve the former unweighted heat-defect majorant of order N⁻³ᐟ²√(Y₀𝒜_N)exp(𝒜_N/4) to

\[
\boxed{\|z_0\|_{\mathcal X(0,\infty)}
\le\frac3{\sqrt2}\sqrt{S_bY_0}\,e^{\mathcal A_N/4}.}          \tag{15}
\]

Thus the extra √𝒜_N loss can be removed. The full sharper expression (13), rather than (14), vanishes at A=0. Amplification exp(𝒜_N/4) remains, so (15) is not an unconditional cutoff-decaying estimate after that dependence is suppressed. This improvement alone also does not replace the accumulated weighted residual in the nonlinear critical error gate: that gate weights the forcing differently.

## 6. Late-time decay from energy alone

Let M_N be the number of retained nonzero lattice modes, and let E_τ^kin=||v(τ)||₂². Exact Galerkin energy decay gives

\[
U(s)^2\le M_NE_\tau^{\rm kin}e^{-2(s-\tau)},\quad
R(s)\le\tfrac12U(s)^2.
\]

For N≥2, λ=N²>2, this yields the fully explicit recent-defect estimate

\[
\boxed{C_{\tau,T}
\le\frac{M_NE_\tau^{\rm kin}}2
\frac{e^{-2\Delta}-e^{-\lambda\Delta}}{\lambda-2}.}             \tag{16}
\]

It proves actual late-time decay for the forcing defect, independently of the earlier accumulated 𝒜_N. An old defect term is still retained through (3) or (4). For a different exact spectral gap λ=2 the quotient is Δe^(−2Δ); for λ<2 the displayed quotient remains nonnegative with its negative denominator. The restriction N≥2 is used only to present a positive denominator in (16).

The coefficient M_N/(N²−2) has size O(N). Consequently the energy bound alone still does not give decay as N→∞ at a fixed time. Choosing late times depending on N does not certify the preceding full PDE comparison.

## 7. A local comparison gate with honest endpoint carry

The high-frequency decay in (2) belongs to the residual defect. The full comparison error can create low frequencies through its interactions with v and with itself. It cannot be assigned the gap N².

Suppose a preceding rigorous comparison has established

\[
\mathfrak h(e(\tau))\le d.
\]

If the reference is changed at τ, d must include the full jump. Let I=[τ,T], V=||v||_{𝒳(I)}, and choose any valid forcing bound

\[
D_R\ge\|z_\tau\|_{\mathcal X(I)},
\quad
D=d+D_R.
\]

For example take √𝓜_(τ,T) from (5), or the minimum of its coarser bounds (6), (7), and (13)/√2; a directly validated signed prefix defect norm can be smaller. If

\[
\gamma=1-2cV>0,\qquad4cD<\gamma^2,
\]

then the existing bilinear heat-space estimate contracts the error ball with radius

\[
\rho=\frac{\gamma-\sqrt{\gamma^2-4cD}}{2c},\qquad
\|e\|_{\mathcal X(I)}\le\rho.                                 \tag{17}
\]

This is a finite-interval certificate, not a reset of e(τ). The initial homogeneous term is bounded by d because the heat norm is defined over the entire future.

A separate endpoint estimate retains decay and the full forcing history. Set

\[
V_T=\left(\int_\tau^T e^{-2(T-s)}U(s)^2ds\right)^{1/2}\le V.
\]

The full mean-zero error uses the true lowest spectral gap 1. Modewise future heat-kernel bounds applied to the exact error equation yield

\[
\boxed{\mathfrak h(e(T))
\le e^{-(T-\tau)}d+\frac{C_{\tau,T}}{\sqrt2}
       +2cV_T\rho+c\rho^2.}                                  \tag{18}
\]

The four terms are respectively the old full error, recent full residual forcing, reference–error interaction, and error self-interaction. For the mixed term the endpoint bilinear estimate is

\[
c\int_\tau^T e^{-(T-s)}U(s)A_0(e(s))ds\le cV_T\rho,
\]

and the factor 2 in the error equation gives 2cV_Tρ. The self term is bounded by c∫exp(−(T−s))A₀(e)²≤cρ². The forcing term alone uses the stronger spectral gap N². Formula (18) is the endpoint bound carried into the next interval. A prefix solution-space radius ρ alone does not automatically bound the endpoint heat norm when D_R was certified only on the prefix. One may additionally use the endpoint bound ρ only if the chosen D_R also bounds the endpoint forcing heat norm; the full-future bound (5) has that property, whereas an arbitrary signed prefix bound need not. No part of an old error is discarded.

The statement can be proved uniformly against every larger Galerkin cutoff before taking a PDE limit: projection contracts all positive Fourier bounds. Every interval preceding T must pass its own gate. If infinitely many intervals are used, their endpoints still require a nonaccumulation argument; short residual memory does not prove it automatically.

## 8. What this improves and what it does not prove

The previous total-budget estimate ignores forcing time. The new results establish three concrete improvements:

1. Old high-output forcing contributions decay on time scale N⁻², with exact endpoint and future-tail estimates (2)–(4).
2. A positive quadratic time-memory budget bounds the entire heat response and gives the factor 1/N for an L²-time residual bound, as in (5)–(6).
3. The covariance identity gives the causal estimate (11) and removes the avoidable √𝒜_N factor from the unweighted global defect bound (15).

For comparison, a nonnegative forcing budget of height r over a duration L has (7) upper bound rL/√2 and (6) upper bound r√L/N. Its exact quadratic majorant is

\[
\mathcal M=r^2\left[{L\over N^2}
              -{1-e^{-N^2L}\over N^4}\right].
\]

Thus (5) interpolates between these small-duration and long-duration bounds while retaining the actual time scale N⁻². Contributions from separated forcing intervals have exponentially damped cross terms in (5). This is an illustrative forcing profile, not a claim that an arbitrary such profile is a Navier–Stokes residual.

These bounds can make a specified finite comparison pass where a whole-time absolute-sum budget is too large. Their hypotheses involve actual reference envelopes, residual moments, and preceding error certificates. No theorem that every smooth datum has such passing intervals or a passing cutoff is obtained. The known unconditional bound 𝒜_N=O(E₀N) remains insufficient for the critical logarithmic bridge. All original NS-B scope restrictions and the absence of a universal regularity proof remain in force.
