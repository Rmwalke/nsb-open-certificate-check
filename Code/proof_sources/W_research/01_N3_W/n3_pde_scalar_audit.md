# Independent PDE and scalar audit for the N=3 certificate

This audit concerns the unforced three-dimensional incompressible Navier–Stokes equations with viscosity 1 on the mean-zero, divergence-free subspace of the 2π-periodic torus. All spatial integrals use the normalized measure `(2π)^−3 dx`. The fixed initial datum is

\[
W(x,y,z)=(\sin z,\sin x,\sin y).
\]

The conclusions below are specific to this datum and to the explicit H¹ neighborhood stated at the end. They make no assertion for arbitrary initial data.

## Computational premises and checks

The finite certificate uses the exact spherical Fourier cutoff N=3, represented in the supplied ten-dimensional invariant basis. Its recorded validation covers [0,8], with dyadic precision 224, step 1/128, and Taylor order 16, followed by an analytic half-line tail estimate. This audit did not run another time integration.

The analytical argument requires the exact Fourier closure and residual identities represented by `n3_exact_algebra_intermediate.json`, together with the mean-zero vector/tensor Sobolev inequality

\[
\|a\|_6\le C_S\|\nabla a\|_2,\qquad C_S=25/3.
\]

This audit independently checked from the authorized exact data: basis reality, mean zero, divergence freedom, orthogonality and the listed squared norms, Laplacian shell values, cubic cancellation in the finite energy identity, the 108 retained wavevectors, and the identity `sum |k|² = 660` over that support. The full Fourier convolution/closure calculation is a separate algebraic premise; these static checks do not replace that calculation. The Sobolev constant is the previously stipulated valid constant, not a claim about the optimal constant.

The reviewed enclosure algorithm uses strict Picard inclusion, outward dyadic arithmetic, rigorous derivative remainders, and exponential variation of constants for endpoint propagation. The residual and gradient bounds are applied to continuous step enclosures. The numerical conclusions used here are the upper bounds

\[
\eta<\frac{587}{10^9},\qquad I<156400,\qquad
\kappa=27C_S^6\eta^2I<\frac{49}{100},\qquad
\Lambda_\infty<\frac{176}{25}.
\]

The stored η and I totals agree exactly with the last recorded finite-interval upper bounds plus the recorded tail upper bounds. The rounded gate is independently checked by rational arithmetic:

\[
27(25/3)^6(587/10^9)^2\,156400
=\frac{134726479}{276480000}<\frac{49}{100}.
\]

The exponential comparison `exp(176/25)<1150` can be checked without floating point. Set a=11/25 and

\[
U=\sum_{j=0}^{16}\frac{a^j}{j!}
+\frac{a^{17}/17!}{1-a/18}.
\]

The positive Taylor series and a geometric tail bound give `exp(a)<U`; exact rational arithmetic gives `U^16<1150`.

## H¹ error inequality

First retain an arbitrary viscosity ν>0 to check the constants. Let v be a smooth divergence-free reference field with residual r, and let e=u−v. Write

\[
z=\|\nabla e\|_2^2,\qquad D=\|\Delta e\|_2^2,
\qquad G=\|\nabla v\|_\infty,
\]

where the pointwise gradient norm is Frobenius. Testing the error equation with −Δe is legitimate for smooth solutions and for finite Galerkin approximations. Leray and Galerkin projections disappear against the test field. Periodicity and divergence freedom give

\[
\frac12 z'+\nu D
\le Gz+G\sqrt z\sqrt D
+C_S^{3/2}z^{3/4}D^{3/4}+\|r\|_2\sqrt D.
\]

Here the first term follows by integrating the transport term `(v·∇)e` once by parts. The second uses `||(e·∇)v||₂≤G||e||₂` and mean-zero Poincaré. The cubic error term follows from

\[
\|e\|_6\|\nabla e\|_3\|\Delta e\|_2
\le C_S^{3/2}z^{3/4}D^{3/4},
\]

using interpolation, the tensor Sobolev inequality for ∇e, and the periodic identity `||D²e||₂=||Δe||₂`.

Allocate νD/4 to each of the latter three terms. In particular,

\[
C_S^{3/2}z^{3/4}D^{3/4}
\le\frac\nu4D+\frac{27C_S^6}{4\nu^3}z^3.
\]

After multiplying by two, this gives exactly

\[
z'+\frac\nu2D
\le\left(2G+\frac{2G^2}{\nu}\right)z
+\frac{27C_S^6}{2\nu^3}z^3+\frac2\nu\|r\|_2^2.
\]

No coefficient in the stated inequality needs enlargement. On this torus `D≥z`.

## Scalar implication of the passing gate

Now set ν=1 and let v=v₃ be the exact N=3 Galerkin reference solution. Define

\[
\Lambda(t)=\int_0^t(2G+2G^2)\,ds,\quad
L(t)=\Lambda(t)-t/2,\quad f(t)=2\|r(t)\|_2^2,
\]

\[
\eta=d^2+\int_0^\infty e^{-L(t)}f(t)\,dt,
\qquad I=\int_0^\infty e^{2L(t)}\,dt,
\]

where d² is the initial squared H¹ error. Setting `w=e^(−L)z` gives

\[
w'\le A e^{2L}w^3+e^{-L}f,
\qquad A=27C_S^6/2.
\]

For completeness, add the remaining forcing budget:

\[
V(t)=w(t)+\int_t^\infty e^{-L(s)}f(s)\,ds.
\]

Then `V(0)=η`, `V≥w≥0`, and `V'≤Ae^(2L)V³`. Bernoulli comparison therefore gives

\[
z(t)\le
\frac{e^{L(t)}\eta}
{\sqrt{1-27C_S^6\eta^2\int_0^t e^{2L(s)}\,ds}}
\le\frac{e^{\Lambda_\infty}\eta}{\sqrt{1-\kappa}}e^{-t/2}.
\]

For the fixed datum W, d=0. The certified bounds imply

\[
z(t)<\frac{1150(587/10^9)}{\sqrt{51/100}}e^{-t/2}
<\frac1{1000}e^{-t/2}.
\]

The last comparison is exact: the square of its coefficient is

\[
\frac{182277001}{204000000000000}<\frac1{10^6}.
\]

The quantity z is the **squared** H¹ seminorm of the error. Accordingly, the corresponding norm bound is `||∇(u−v₃)||₂ < (1/√1000)e^(−t/4)`.

## The analytic tail and the saved Λ bound

Let E₈ be the validated upper bound for `||v₃(8)||₂²`, and write τ=t−8. The exact Galerkin energy law and Poincaré give

\[
\|v_3(8+\tau)\|_2^2\le E_8e^{-2\tau}.
\]

The support constants checked above imply

\[
G(8+\tau)\le\sqrt{660E_8}\,e^{-\tau},\qquad
\|r(8+\tau)\|_2^2\le972E_8^2e^{-4\tau}.
\]

The second estimate uses `||v₃||∞≤√108||v₃||₂`, `||∇v₃||₂≤3||v₃||₂`, and contraction of the orthogonal residual projection. Hence

\[
\Lambda_\infty\le\Lambda(8)+2\sqrt{660E_8}+660E_8,
\]

\[
\eta_{\rm tail}\le\frac47\,972E_8^2e^{-L(8)},\qquad
I_{\rm tail}\le e^{2L(8)+2(2\sqrt{660E_8}+660E_8)}.
\]

The saved Λ∞ upper bound correctly uses the sharper expression
`Λ(8)_upper + 2 ceil_Q(sqrt(660 E₈)) + 660 E₈`.
Its square term is exact; it does not need to square the rounded square root. This explains the one-dyadic-unit difference from the slightly less sharp interval-square calculation. Both calculations give `Λ∞<176/25`.

## Uniform higher-Galerkin estimates and the PDE limit

The argument need not assume beforehand that a global PDE solution exists. For every spherical cutoff M≥3, let u_M be the exact Galerkin solution with initial datum W. Finite Galerkin solutions exist globally by their energy identity. Since v₃ belongs to the M-mode space, the difference `e_M=u_M−v₃` satisfies the preceding estimate with residual `P_M r`; its norm is bounded by `||r||₂`. The initial error is zero. Thus the same scalar certificate gives, uniformly in M,

\[
\|\nabla e_M(t)\|_2^2<10^{-3}e^{-t/2}.
\]

The reference solution has spectral support `|k|≤3`, and its energy obeys

\[
\|v_3(t)\|_2^2\le\frac32e^{-2t},\qquad
\|\nabla v_3(t)\|_2^2\le\frac{27}{2}e^{-2t}.
\]

Consequently, with `X_M=||∇u_M||₂²`,

\[
\sqrt{X_M(t)}
<\sqrt{27/2}\,e^{-t}+\sqrt{1/1000}\,e^{-t/4},
\qquad \sup_{M\ge3,t\ge0}X_M(t)<14.
\]

The bound follows, for example, by squaring the right side at t=0 and checking `(√(27/2)+√(1/1000))²<14`. The exact Galerkin energy law gives

\[
\int_0^\infty X_M(t)\,dt\le\frac34,
\qquad
\int_0^\infty X_M(t)^2\,dt<\frac{21}{2}.
\]

These statements apply to all higher Galerkin approximations, not merely to the ten-dimensional reference trajectory.

For cutoffs M<3, the direct spectral energy bound gives `X_M(t)≤9||P_M W||₂² e^(−2t)≤27/2<14`; the energy integral remains at most 3/4. Thus the same displayed supremum and integral bounds hold uniformly over all cutoffs. The uniform constant for M≥3 is itself strictly below 14, so taking the supremum over M preserves the strict margin.

There is also a uniform H² estimate. Integrate the original error inequality, retaining its positive `D/2` term. All right-hand integrals are bounded independently of M: the validated tail controls G and f, while the displayed error bound controls z and z³. Thus `e_M` is uniformly bounded in `L²(0,T;H²)` for each finite T (indeed on the whole half-line), and `u_M` is uniformly bounded in `L∞(0,T;H¹) ∩ L²(0,T;H²)`.

The equation then bounds `∂_t u_M` in `L²(0,T;L²)`. One direct estimate is

\[
\|(u_M\cdot\nabla)u_M\|_2^2
\le C_S^3X_M^{3/2}\|\Delta u_M\|_2
\le\frac{C_S^3}{2}(X_M^3+\|\Delta u_M\|_2^2).
\]

Standard compactness therefore produces a global strong solution in `C([0,T];H¹) ∩ L²(0,T;H²)` for each finite T; strong local convergence passes the nonlinear term. The initial datum is attained, and weak lower semicontinuity followed by time continuity transfers the error bound to the limit. Strong-solution uniqueness makes the limit independent of the selected subsequence. Equivalently, the same uniform H¹ estimate prevents finite-time breakdown of the maximal local strong solution. Since W is smooth, parabolic regularity and continuation give a global smooth solution.

This is a continuation/compactness argument for the fixed datum, not an assumption of general three-dimensional regularity.

## Explicit H¹ neighborhood

Let h be mean-zero and divergence-free, with

\[
\|\nabla h\|_2\le\frac1{2500}.
\]

Keep the same reference solution v₃ and its same residual and growth bounds. For initial datum `W+h`, the initial error is at most `(1/2500)²`; in the M-mode argument use `P_M h`, whose H¹ seminorm is no larger. Therefore

\[
\eta_h<\frac{587}{10^9}+\frac1{2500^2}
=\frac{747}{10^9},
\]

\[
\kappa_h<27(25/3)^6(747/10^9)^2\,156400
=\frac{8080797}{10240000}<\frac45.
\]

The scalar estimate yields the convenient uniform bound

\[
\|\nabla(u-v_3)(t)\|_2^2<\frac1{500}e^{-t/2},
\]

because `[1150(747/10^9)]²/(1−4/5)<(1/500)²`. The preceding strong-solution construction and continuation argument apply. Thus every initial datum in this explicit closed H¹ ball has a unique global strong solution, smooth for positive time; smooth perturbations give global smooth solutions with their smooth initial regularity.

The constants `sup X_M<14` and `∫X_M²<21/2` established above were stated for the fixed datum W. The perturbation-ball conclusion does not rely on silently carrying over its unchanged initial energy.

## Scope and remaining proof obligations

No gap was found in the H¹ testing constants, scalar comparison, tail formulas, requested rational margins, or passage from the uniform Galerkin estimate to the PDE conclusion. These implications use the exact finite-model/residual identities and the valid Sobolev constant as explicit premises. This written audit is not a machine-checked proof and is not a second execution of the numerical enclosure run; reproducibility requires retaining the exact algebra, interval-core and runner versions, and the recorded upper bounds together.

The result is the fixed-W theorem and the explicit radius-1/2500 H¹ neighborhood theorem. Regularity for arbitrary three-dimensional initial data remains open. There is no novelty claim here.
