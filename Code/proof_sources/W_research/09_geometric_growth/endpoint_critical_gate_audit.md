# Independent audit of the endpoint critical reference gate

This new note checks the proposed endpoint transformation and completed-square residual bound. It uses the normalized three-dimensional periodic Navier--Stokes setting, viscosity one, real mean-zero divergence-free fields, and the improved pair estimate in `pair_geometry_adversarial_audit.md`. No historical script, numerical time integration, web search, upload, or modification of prior notes was used.

**Finding:** the endpoint construction is valid. It replaces the degenerating constant-barrier choice at exponent \(1/4\) by the increasing function \(\Psi(b)=be^{-b/2}\), and it gives a finite residual estimate with a polynomial factor in the actual reference amplification budget. The final gate remains sufficient and conditional.

## 1. Critical error inequality before constant Young splitting

Let \(v\) be a fixed smooth finite Galerkin reference, and let \(e=u-v\). Write
\[
b=A_{-1}(e),\qquad h=A_1(e),\qquad U=A_0(v),
\]
\[
R=A_{-1}(v'+Av+B(v,v)),\quad
\mathcal A(t)=\int_0^tU(s)^2\,ds,\quad
P=\mathcal A(\infty)<\infty.
\]
The critical coefficient estimate gives the unsplit inequality
\[
b'+(1-b/2)h\le U\sqrt{bh}+R.                       \tag{1}
\]
Indeed, the self-interaction contributes at most \(A_0(e)^2/2\le bh/2\), while the two cross terms together contribute at most \(U A_0(e)\le U\sqrt{bh}\). For finite Galerkin comparisons, these are inequalities for absolutely continuous norms, interpreted almost everywhere or with upper derivatives at zero coefficients.

While \(0\le b<2\), put \(d=1-b/2>0\). Maximizing the right side of
\[
b'\le U\sqrt b\sqrt h-dh+R
\]
over \(\sqrt h\ge0\) yields
\[
b'\le\frac{U^2b}{4d}+R.                            \tag{2}
\]

## 2. The monotone transform and first-contact argument

Define
\[
\Psi(b)=be^{-b/2},\qquad0\le b\le2.
\]
Then
\[
\Psi'(b)=e^{-b/2}(1-b/2),
\]
which is positive on \([0,2)\), nonnegative on the closed interval, and at most one there. Multiplication of (2) by this nonnegative derivative gives, with \(z=\Psi(b)\),
\[
z'\le\frac{U^2}{4}z+\Psi'(b)R
\le\frac{U^2}{4}z+R.                               \tag{3}
\]
Thus, as long as \(b<2\),
\[
z(t)\le e^{\mathcal A(t)/4}
\left[\Psi(b_0)+\int_0^t e^{-\mathcal A(s)/4}R(s)\,ds\right].
\]
If
\[
\boxed{b_0<2,\qquad
M:=e^{P/4}\left[\Psi(b_0)+
\int_0^\infty e^{-\mathcal A(t)/4}R(t)\,dt\right]
<\frac2e,}                                         \tag{4}
\]
then no first contact with \(b=2\) is possible: that contact would require \(z=\Psi(2)=2/e\), contradicting the same strict upper bound \(M\). Since \(\Psi\) is strictly increasing before two, its inverse on \([0,2/e)\) gives
\[
\boxed{b(t)\le B:=\Psi^{-1}(M)<2.}                 \tag{5}
\]
The inverse is the branch taking values in \([0,2)\). No global monotonicity beyond that interval is needed. In particular, the argument must first verify \(b_0<2\); small \(\Psi(b_0)\) alone would not imply it for arbitrary large \(b_0\).

## 3. Dissipation and continuation are retained

Put \(\delta=1-B/2>0\). Return to (1), rather than trying to infer dissipation from (3). Since \(b\le B\),
\[
b'+\delta h\le U\sqrt{Bh}+R
\le\frac\delta2h+\frac{B}{2\delta}U^2+R.
\]
Consequently
\[
\boxed{\int_0^\infty h(t)\,dt
\le H_*:=\frac{2b_0}{\delta}
+\frac{BP}{\delta^2}
+\frac2\delta\int_0^\infty R(t)\,dt.}             \tag{6}
\]
For an exact finite Galerkin reference, its full exterior residual satisfies \(R\le U^2/2\), so \(\int R\le P/2\). Thus all terms on the right are finite and explicit from the reference quantities. The endpoint transform has not sacrificed the dissipation needed for continuation.

In particular
\[
\int_0^\infty A_0(e)^2\,dt\le BH_*,\qquad
\|u\|_{L^2_tA_0}\le\sqrt P+\sqrt{BH_*}.           \tag{7}
\]
Applying the comparison uniformly to all larger Galerkin cutoffs and then the earlier \(H^1\) estimate proves global strong continuation and, for smooth initial data, global smoothness. It also gives
\[
\sup_{L\ge N}\int_0^\infty\|\nabla u_L(t)\|_2^4dt
\le\frac{E_0X_0}{2}
\exp\left[(\sqrt P+\sqrt{BH_*})^2\right].          \tag{8}
\]
The finitely many cutoffs \(L<N\) may be included by taking the maximum with \(N^2E_0^2/2\), using their energy identity and \(\|\nabla u_L\|_2^2\le N^2E_0\). Thus the endpoint criterion supplies an actual NS-B bound, not only an upper bound on a transformed error.

## 4. The endpoint residual estimate

For the spherical reference \(v_N\), write
\[
Y=\|\Lambda^2v_N\|_2^2,\qquad
Z=\|\Lambda^3v_N\|_2^2,
\qquad Y_0=\|\Lambda^2u_0\|_2^2.
\]
The improved coefficient-one derivative estimate gives
\[
Y'+2Z\le2U\sqrt{YZ}.                               \tag{9}
\]
Set
\[
E=e^{-\mathcal A/2}Y,\qquad W=e^{-\mathcal A/2}Z.
\]
Direct differentiation, including the derivative of the weight, proves
\[
\boxed{E'+2\left(\sqrt W-\tfrac12U\sqrt E\right)^2\le0.} \tag{10}
\]
It follows that \(E(t)\le Y_0\) and
\[
\int_0^\infty
\left(\sqrt W-\tfrac12U\sqrt E\right)^2dt\le Y_0/2.
\tag{11}
\]
The exterior input-band localization is
\[
R_N\le\sqrt{S_6^{\rm band}(N)}\,U\sqrt Z,
\qquad S_6^{\rm band}(N)=\sum_{N/2<|k|\le N}|k|^{-6}.
\]
After multiplying by the residual weight \(e^{-\mathcal A/4}\), this becomes \(\sqrt{S_6^{\rm band}}\,U\sqrt W\). Decompose
\[
\sqrt W=\tfrac12U\sqrt E+
\left(\sqrt W-\tfrac12U\sqrt E\right).
\]
The term in parentheses may have either sign. For an upper bound its integral against \(U\) is bounded by Cauchy--Schwarz using (11), with no pointwise positivity assumed. Therefore
\[
\begin{aligned}
\int_0^\infty e^{-\mathcal A/4}R_N\,dt
&\le\sqrt{S_6^{\rm band}(N)}
\left[\tfrac12\sqrt{Y_0}\int_0^\infty U^2dt
+\left(\int_0^\infty U^2dt\right)^{1/2}
 \left(\frac{Y_0}{2}\right)^{1/2}\right]\\
&\le\boxed{\sqrt{S_6^{\rm band}(N)Y_0}
\left(\frac P2+\sqrt{\frac P2}\right).}
\end{aligned}
\tag{12}
\]
This estimate requires the actual finite value \(P\); it is not uniform in an unknown unbounded reference amplification budget. It is nevertheless finite at the endpoint exponent \(1/4\), unlike the earlier denominator bound as \(\beta\downarrow1/4\).

## 5. A direct gate using full-data norms

The initial truncation error obeys
\[
b_0\le\sqrt{Y_0S_6^{\rm tail}(N)},\qquad
S_6^{\rm tail}(N)=\sum_{|k|>N}|k|^{-6}.
\]
Also \(\Psi(b_0)\le b_0\) for nonnegative \(b_0\). Define the conservative majorant
\[
\overline M_N=
e^{P/4}\sqrt{Y_0}\left[
\sqrt{S_6^{\rm tail}(N)}+
\left(\frac P2+\sqrt{\frac P2}\right)
\sqrt{S_6^{\rm band}(N)}\right].
\]
Then the fully explicit sufficient test is
\[
\boxed{\overline M_N<2/e.}                          \tag{13}
\]
This condition also verifies the separate initial-side condition in (4), since
\[
b_0\le\sqrt{Y_0S_6^{\rm tail}(N)}
\le e^{-P/4}\overline M_N<2/e<2.
\]
Once that condition is established, the transform argument applies and its \(M\) is at most \(\overline M_N\). A certified actual \(b_0\) and actual weighted residual may of course pass (4) even when the simpler majorant (13) fails.

The lattice factors in (13) are of order \(N^{-3/2}\); its extra factor is of order \(1+P\). Hence a strict asymptotic slope below six still makes this endpoint test pass along the corresponding unbounded subsequence. Nothing here proves that such amplification control holds for arbitrary smooth data. Equality cases of (4) or (13) are not included.

**Audit status:** no fatal gap was found in the endpoint transform, first-contact barrier, retained dissipation, completed-square residual bound, or the combined finite gate. The exponent \(1/4\) is reached by a different valid barrier and a finite polynomial residual bound, not by substituting an inadmissible parameter value into a divergent earlier formula.
