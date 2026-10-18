# A heat-budget fixed point in the Fourier-amplitude space

This is a written estimate, with no new ODE/PDE integration. It concerns the unforced, mean-zero, real, divergence-free Navier–Stokes equation with viscosity one on the normalized \(2\pi\)-torus. Every Fourier sum below uses all nonzero lattice vectors and the Euclidean norm of the complex vector coefficient. No novelty claim is made.

## 1. Space, heat norm, and the Leray factor

Write
\[
 A_s(w)=\sum_{k\ne0}|k|^s|\widehat w_k|,
 \qquad
 \|w\|_{\mathcal X(I)}=
 \left(\int_I A_0(w(t))^2\,dt\right)^{1/2},
 \qquad \mathcal X=\mathcal X(0,\infty).
\]
This is the Banach space \(L^2_t\ell^1_k\), restricted to the closed real, divergence-free, mean-zero subspace. Let \(S(t)=e^{t\Delta}\), \(A=-\Delta\), and define the symmetric bilinear nonlinearity
\[
 B(u,v)=\tfrac12\mathbb P\operatorname{div}
       (u\otimes v+v\otimes u),\qquad
 Q(u,v)(t)=\int_0^t S(t-s)B(u,v)(s)\,ds.
\]
Thus the mild equation is \(u=S(t)u_0-Q(u,u)\).

For a real unit vector \(n\), let \(a=\alpha n+a_\perp\) and \(b=\beta n+b_\perp\), with complex coefficients. Then
\[
 |\mathbb P_{n^\perp}[(n\cdot a)b+(n\cdot b)a]|
 =|\alpha b_\perp+\beta a_\perp|
 \le |\alpha||b_\perp|+|\beta||a_\perp|
 \le |a||b|.
\]
The last step is the two-dimensional Cauchy inequality for the two nonnegative pairs. Consequently, after reindexing the second term in the symmetric Fourier convolution,
\[
 |\widehat {B(u,v)}_k|
 \le\frac{|k|}{2}\sum_{p+q=k}|\widehat u_p||\widehat v_q|,
 \qquad k\ne0.                                      \tag{1}
\]
The zero output mode is exactly zero because of the derivative. It must be omitted; no zero-frequency heat-kernel inverse is used. Formula (1) is valid for mixed \(u,v\), not only on the diagonal. The resulting operator constant below is a proved upper bound, not a claim of optimality for the full operator.

## 2. The global bilinear estimate is valid

For each nonzero \(k\), put \(K_k(t)=\mathbf1_{t\ge0}|k|e^{-|k|^2t}\). Its exact norm is
\[
 \|K_k\|_{L^2(0,\infty)}^2
 =\int_0^\infty |k|^2e^{-2|k|^2t}\,dt=\frac12.
\]
Apply Minkowski to the Fourier sum and then Young's \(L^1_t\to L^2_t\) convolution estimate, keeping this individual kernel until that step. Tonelli applies to the nonnegative majorants. This gives
\[
\begin{aligned}
 \|Q(u,v)\|_{\mathcal X}
 &\le\frac1{2\sqrt2}
       \sum_{k\ne0}\sum_{p+q=k}
          \int_0^\infty |\widehat u_p(s)||\widehat v_q(s)|\,ds\\
 &\le\frac1{2\sqrt2}
       \int_0^\infty A_0(u(s))A_0(v(s))\,ds\\
 &\le c\|u\|_{\mathcal X}\|v\|_{\mathcal X},
 \qquad \boxed{c=\frac1{2\sqrt2}}.                    \tag{2}
\end{aligned}
\]
All exchanges can first be made for finite Fourier sums and compactly supported time functions, then extended by density and completeness. The same constant works on every finite interval, for time-shifted operators, and with an orthogonal spherical Fourier projection on the output. Replacing the modewise kernel by its pointwise supremum before integration would lose this strong endpoint estimate; that replacement is unnecessary.

For a datum \(a\), define its heat norm by
\[
 \mathfrak h(a)=\|S(t)a\|_{\mathcal X},\qquad
 P_{\rm heat}(a)=\mathfrak h(a)^2
 =\sum_{k,\ell\ne0}
   \frac{|\widehat a_k||\widehat a_\ell|}
        {|k|^2+|\ell|^2}.                            \tag{3}
\]
The double sum is an exact identity by Tonelli, not a numerical approximation. It is a norm after taking the square root. The heat semigroup and spherical Fourier projections contract this norm. The familiar estimate
\[
 P_{\rm heat}(a)\le\tfrac12 A_{-1}(a)^2             \tag{4}
\]
follows from \(|k|^2+|\ell|^2\ge2|k||\ell|\), but is often much less precise than (3).

## 3. Small heat budget gives a global solution

Set \(H=\mathfrak h(u_0)\). If
\[
 \boxed{P_{\rm heat}(u_0)<\tfrac12},                 \tag{5}
\]
then \(4cH<1\). Define
\[
 r_-=\frac{1-\sqrt{1-4cH}}{2c}
     =\sqrt2\bigl(1-\sqrt{1-\sqrt2H}\bigr).
                                                               \tag{6}
\]
The mild map sends the closed radius-\(r_-\) ball in \(\mathcal X\) into itself, since \(H+cr_-^2=r_-\), and has contraction factor \(2cr_-<1\). Hence there is a global mild solution with
\[
 \|u\|_{\mathcal X}\le r_-\le2H,\qquad r_-^2<2.     \tag{7}
\]
The endpoint equality \(P_{\rm heat}=1/2\) is not covered by this argument. Uniqueness also holds among arbitrary solutions with finite \(\mathcal X\) norm, by subdividing time into intervals on which \(c(\|u\|_{\mathcal X(I)}+\|v\|_{\mathcal X(I)})<1\), carrying the identical past forward. Finite total squared norms allow a finite such subdivision, with a final infinite interval.

For smooth data, the conclusion is global smooth regularity, not merely a mild distributional construction. Here is a direct uniform Galerkin proof. Every finite spherical Galerkin solution \(u_N\) exists globally by its exact energy equality. Its prefix norm
\[
 y_N(T)=\|u_N\|_{\mathcal X(0,T)}
\]
is continuous, starts at zero, and (2) gives \(y_N(T)\le H+cy_N(T)^2\). The two roots of the quadratic are separated under (5). Continuity keeps \(y_N(T)\) in the lower branch, so \(y_N(T)\le r_-\) for every \(N,T\).

Put \(X_N=\|\nabla u_N\|_2^2\), \(Y_N=\|\Delta u_N\|_2^2\), \(U_N=A_0(u_N)\). Since \(\|u_N\|_\infty\le U_N\), direct testing by \(-\Delta u_N\) gives
\[
 \tfrac12X_N'+Y_N
 \le U_N\sqrt{X_NY_N}
 \le\tfrac12Y_N+\tfrac12U_N^2X_N,
 \quad\text{hence}\quad X_N'+Y_N\le U_N^2X_N.       \tag{8}
\]
An integrating factor and (7) yield the explicit uniform bounds
\[
 \boxed{\sup_{t\ge0}X_N(t)\le X_0e^{r_-^2},\qquad
        \int_0^\infty Y_N(t)\,dt\le X_0e^{r_-^2}},  \tag{9}
\]
where \(X_0=\|\nabla u_0\|_2^2\). In particular both right-hand sides are at most \(X_0e^{4P_{\rm heat}(u_0)}\le e^2X_0\). These estimates give global strong compactness and continuation in \(H^1\); the usual differentiated local estimates give smoothness for smooth data. Thus no global regularity premise and no external Serrin criterion are needed. For \(H^1\) data, the same proof gives a global strong solution, smooth at positive times.

There is also an explicit bound for the requested squared-enstrophy integral. With \(E_0=\|u_0\|_2^2\), the exact Galerkin energy equality gives \(\int_0^\infty X_N\le E_0/2\). Combining this with (9),
\[
 \boxed{\int_0^\infty X_N(t)^2\,dt
 \le \frac{E_0X_0}{2}e^{r_-^2}
 \le \frac{E_0X_0}{2}e^{4P_{\rm heat}(u_0)}}.          \tag{9a}
\]
For nonzero data this is strictly less than \(4E_0X_0\), using \(4P_{\rm heat}<2\) and \(e^2<8\); for zero data every quantity vanishes. The same bounds hold for the limiting solution. These are cutoff-independent estimates for this explicit data class. They use no Sobolev constant \(C_S\) and no N3 numerical certificate.

For any perturbation family independently certified to satisfy \(P_{\rm heat}(u_0)<49/100\), one may use
\[
 r_-^2<
 2\left(1-\sqrt{1-\frac{7\sqrt2}{10}}\right)^2
 \le \frac{49}{25},
\]
or simply the upper exponent \(49/25\) in (9) and (9a), with each datum's \(E_0,X_0\).

## 4. The criterion need not be small in other quoted norms

Neither a small \(H^{1/2}\) norm nor a fixed upper bound on \(A_{-1}\) is necessary for (5). A concrete smooth finite family makes this explicit. Let \(q_j\) cycle through the three coordinate unit vectors, and let \(e_j\) be the next coordinate unit vector, so \(e_j\cdot q_j=0\). Define
\[
 u_0^{(J)}(x)=\frac{\varepsilon}{\sqrt J}
       \sum_{j=0}^{J-1}2^j e_j\sin(2^jq_j\cdot x).
\]
The total Fourier amplitude on shell \(2^j\) is \(\varepsilon2^j/\sqrt J\), so
\[
 A_{-1}(u_0^{(J)})=\varepsilon\sqrt J,
\quad
 P_{\rm heat}(u_0^{(J)})
 \le\varepsilon^2\left[\frac12+
              2\sum_{d=1}^\infty\frac1{2^d+2^{-d}}\right]
 <\frac{23}{10}\varepsilon^2.                       \tag{10}
\]
For the last bound, the \(d=1\) term is \(2/5\), and the sum over \(d\ge2\) is strictly less than \(\sum_{d\ge2}2^{-d}=1/2\). Taking \(\varepsilon=1/4\) gives \(P_{\rm heat}<23/160<1/2\), uniformly in \(J\), while \(A_{-1}\to\infty\). Also
\[
 \|u_0^{(J)}\|_{\dot H^{1/2}}^2
 =\frac{\varepsilon^2}{2J}\sum_{j=0}^{J-1}2^{3j}\longrightarrow\infty.
\]
For \(J\ge3\) these fields involve all three coordinates and have nonzero nonlinearity: the first two modes alone generate the nonzero \((1,2,0)\) output in the third component. This is a statement about the scope of the norm criterion, not a novelty claim or a claim that all large data pass it. A small Sobolev norm can still be a sufficient way to bound the heat budget from above. In the opposite direction, the previously certified datum \(W=(\sin z,\sin x,\sin y)\) has \(P_{\rm heat}(W)=9/2\) and fails (5), despite its separate N3 certificate. This gate is sufficient and not necessary.

## 5. A reference-based gate retains the full residual

Let \(v(t)\) be a smooth, real, divergence-free, mean-zero finite Fourier reference, with residual
\[
 r=v'+Av+B(v,v),\qquad e_0=u_0-v(0).
\]
The exact error equation in mild form is
\[
 e=g-2Q(v,e)-Q(e,e),\qquad
 g(t)=S(t)e_0-\int_0^tS(t-s)r(s)\,ds.               \tag{11}
\]
There is no additional low-low forcing: it is already part of the full residual \(r\). Put \(V=\|v\|_{\mathcal X}\), \(D=\|g\|_{\mathcal X}\). A sufficient all-time gate is
\[
 \boxed{2cV<1,\qquad 4cD<(1-2cV)^2}.              \tag{12}
\]
Indeed the error map contracts the ball of radius
\[
 \rho=\frac{1-2cV-
       \sqrt{(1-2cV)^2-4cD}}{2c},                 \tag{13}
\]
with Lipschitz constant \(2cV+2c\rho<1\). An independently validated norm of the signed defect \(g\) can retain cancellation and heat damping. A simpler entirely positive bound is
\[
 D\le\mathfrak h(e_0)+\frac1{\sqrt2}
                         \int_0^\infty A_{-1}(r(s))\,ds,        \tag{14}
\]
because \(\|\mathbf1_{t\ge s}e^{-|k|^2(t-s)}\|_{L^2_t}=1/(\sqrt2|k|)\). The zero residual mode must vanish; that is automatic for a mean-zero reference in this setting.

If \(u_0\in H^1\), (12) yields global strong regularity and the explicit analogues of (9) and (9a) with \(e^{(V+\rho)^2}\) instead of \(e^{r_-^2}\). To see the uniform approximation directly, for each higher Galerkin cutoff containing the reference support the signed defect is the projection of \(g\), which contracts \(\mathcal X\). The continuous prefix-norm quadratic barrier gives the same error bound. Hence (8) applies uniformly to the full solutions with \(\int A_0(u_N)^2\le(V+\rho)^2\).

## 6. Large finite reference norm: finite interval gates and honest carryover

The restriction \(2cV<1\) in (12) is avoidable by a finite partition, provided the full reference norm is finite. This yields an explicit, conservative stability neighborhood of any such reference.

Let \(0=t_0<t_1<\cdots<t_m=\infty\), with finite intermediate endpoints. On \(I_i=[t_{i-1},t_i]\), set
\[
 V_i=\|v\|_{\mathcal X(I_i)},\qquad
 R_i=\frac1{\sqrt2}\int_{I_i}A_{-1}(r(s))\,ds.
\]
If \(\mathfrak h(e(t_{i-1}))\le d_i\), put \(D_i=d_i+R_i\). Provided
\[
 2cV_i<1,\qquad4cD_i<(1-2cV_i)^2,
\]
formula (13), with \(V_i,D_i\), gives an interval error norm bound \(\rho_i\). At a finite right endpoint,
\[
 \mathfrak h(e(t_i))
 \le d_i+R_i+2cV_i\rho_i+c\rho_i^2=\rho_i.          \tag{15}
\]
For (15), evolve each endpoint Duhamel coefficient into future heat time. Its kernel has norm
\(
 |k|\,\|e^{-|k|^2(t_i-s+\tau)}\|_{L^2_{\tau>0}}
 =e^{-|k|^2(t_i-s)}/\sqrt2\le1/\sqrt2
\).
The residual uses the analogous kernel without \(|k|\). Thus \(d_{i+1}=\rho_i\) is an honest full-support carryover; it does not reset the error. Sharper endpoint enclosures or directly computed signed defects may improve it.

Here is one fully explicit conservative version. Put \(P_v=V^2<\infty\) and choose
\[
 m=\max\{1,\lceil2P_v\rceil\},\qquad V_i^2\le\tfrac12.
\]
A finite partition with these bounds exists by absolute continuity of the cumulative integral; it has no accumulating finite endpoints. If
\[
 F=\mathfrak h(e_0)+\frac1{\sqrt2}\int_0^\infty A_{-1}(r(s))\,ds
 \quad\text{satisfies}\quad
 \boxed{F<\frac{1}{4^{m+1}c}},                     \tag{16}
\]
then every interval gate passes. Indeed \(2cV_i\le1/2\), the threshold \(D_i<1/(16c)\) suffices, and rationalizing (13) gives \(\rho_i\le4D_i\). More precisely, \(D_i\le4^{i-1}(d_1+\sum_{j\le i}R_j)\le4^{i-1}F\), and \(\rho_i\le4^iF\); (16) makes the strict threshold hold even on the last interval. In particular
\[
 \|e\|_{\mathcal X}^2
 \le F^2\sum_{i=1}^m16^i
 =\frac{16(16^m-1)}{15}F^2.                        \tag{17}
\]
This is a usable existence bound for arbitrary finite \(P_v\), not an efficient numerical constant. The individual gates and signed defect norms are preferable in applications. The argument does not claim that multiplication by a decreasing scalar time weight automatically decreases the universal constant in (2); simple weighted Minkowski estimates retain the undiminished local kernel near its starting time.

The same construction works on a finite horizon for any smooth finite reference, without requiring finite all-time \(P_v\). A finite stopping rule is
\[
 \mathfrak h(v(T))+d_{\rm end}<1/\sqrt2,             \tag{18}
\]
where \(d_{\rm end}\) is a certified bound for the complete endpoint error heat norm. Then (5), restarted at \(T\), proves global continuation. Every finite interval before \(T\) must have passed its own gate. Neither selecting a new reference nor increasing its cutoff erases the old error; a jump must be included in the new heat-norm error bound unless a tighter independent recentering enclosure is available.

## Scope

The universal bilinear bound (2), strict heat-budget criterion (5), uniform smooth-data continuation bounds (9), and residual-based stability gates are established by the written estimates above. No arbitrary-data success follows: an arbitrary smooth datum need not satisfy (5), and neither residual smallness nor the interval gates have been proved uniformly for all data. The Fourier-amplitude solution-space argument is a standard type of mild fixed-point method; these calculations make the constants and the full-residual accounting explicit without asserting novelty.
