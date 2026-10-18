# A sufficient convergence criterion for repeated residual heat corrections

Private written NS-B continuation. This note uses no old program or import, evolution computation, quadrature, external search, upload, or modification of earlier files. It establishes a local a posteriori theorem, including a version for references with arbitrarily large finite time-integrated Fourier amplitude. It does not assert that arbitrary data meet its small-defect condition.

Use viscosity one, normalized measure on the 2π torus, and real mean-zero divergence-free fields. Put
\[
\mathcal X(I)=L^2(I;A_0),\quad
\mathcal Q(u,v)(t)=\int_0^t e^{(t-s)\Delta}B_s(u,v)(s)\,ds,
\quad c=\frac1{2\sqrt2},
\]
where \(B_s\) is the symmetric NS bilinear operator. The previously proved estimate is
\[
\|\mathcal Q(u,v)\|_{\mathcal X(0,\infty)}
\le c\|u\|_{\mathcal X}\|v\|_{\mathcal X}.              \tag{1}
\]
It remains valid for zero-extended interval inputs and for projections on the output. Write \(S(t)=e^{t\Delta}\), and \(\mathcal H r(t)=\int_0^tS(t-s)r(s)\,ds\).

## 1. The exact iteration and its full residual

Let \(v_0\) be a smooth finite Fourier reference with initial value \(a\), finite \(\mathcal X\) norm, and full residual
\[
r_0=v_0'-\Delta v_0+B_s(v_0,v_0).
\]
For every stage define
\[
z_j=\mathcal Hr_j,\qquad v_{j+1}=v_j-z_j,\qquad
r_j=v_j'-\Delta v_j+B_s(v_j,v_j).
\]
The identity
\[
\mathcal Hr_j=v_j-Sa+\mathcal Q(v_j,v_j)
\]
shows exactly that
\[
\boxed{v_{j+1}=Sa-\mathcal Q(v_j,v_j).}             \tag{2}
\]
Thus the repeated heat-defect correction is the usual mild Picard map started from a chosen reference. It is not Newton iteration and does not have an automatic quadratic convergence rate.

Every correction preserves \(v_j(0)=a\). Its new residual is exactly
\[
\boxed{r_{j+1}=-2B_s(v_j,z_j)+B_s(z_j,z_j).}       \tag{3}
\]
All low output modes are retained. If \(v_0\) has support in \(|k|\le K\), then each finite iterate has support in \(|k|\le2^jK\), and its residual in \(|k|\le2^{j+1}K\). No high-frequency lower support bound persists in general.

For \(V_j=\|v_j\|_{\mathcal X}\), \(D_j=\|z_j\|_{\mathcal X}\), (1) gives
\[
V_{j+1}\le V_j+D_j,\qquad
D_{j+1}\le c(V_j+V_{j+1})D_j
\le 2cV_jD_j+cD_j^2.                              \tag{4}
\]
The presence of the linear term matters: without a suitable contraction estimate, repeated correction has no guaranteed convergence.

## 2. A simple global-ball test

Set \(V=\|v_0\|_{\mathcal X}\), \(D=\|\mathcal Hr_0\|_{\mathcal X}\). If
\[
\gamma=1-2cV>0,\qquad 4cD<\gamma^2,                \tag{5}
\]
define
\[
\rho=\frac{\gamma-\sqrt{\gamma^2-4cD}}{2c}.
\]
The map \(T(f)=Sa-\mathcal Q(f,f)\) preserves the closed \(\mathcal X\) ball of radius \(\rho\) around \(v_0\), since
\[
\|T(v_0+e)-v_0\|\le D+2cV\|e\|+c\|e\|^2.
\]
It contracts there with
\[
q=2c(V+\rho)<1.
\]
Consequently all iterations (2) converge to the solution with initial value \(a\), and
\[
\|v_{j+1}-v_j\|_{\mathcal X}\le q^jD,\qquad
\|u-v_j\|_{\mathcal X}\le\frac{q^jD}{1-q}.          \tag{6}
\]
This makes explicit the convergence content of the existing heat-space reference gate. It is useful when \(V<\sqrt2\), but it cannot by itself certify a reference with larger \(V\).

## 3. A causal block norm removes the small-reference restriction

Take a finite time partition
\[
0=t_0<t_1<\cdots<t_{m-1}<t_m=\infty,\qquad I_i=[t_{i-1},t_i].
\]
Choose positive weights \(0<w_i\le1\), normalized so their maximum is one. A common rescaling always permits this normalization. For an error \(e\), write
\[
E_i=\|e\|_{\mathcal X(I_i)},\qquad
\|e\|_*=\max_{1\le i\le m}w_iE_i,\qquad
V_i=\|v_0\|_{\mathcal X(I_i)}.
\]
There are finitely many blocks, so this is an equivalent Banach norm:
\[
\|e\|_{\mathcal X}\le K_w\|e\|_*,
\qquad K_w=\left(\sum_{i=1}^m w_i^{-2}\right)^{1/2}.            \tag{7}
\]

Causality and (1) give, for every output block \(i\),
\[
\|\mathcal Q(u,v)\|_{\mathcal X(I_i)}
\le c\sum_{j\le i}\|u\|_{\mathcal X(I_j)}
                       \|v\|_{\mathcal X(I_j)}.               \tag{8}
\]
To prove this, split the input forcing into time blocks. Each block's complete future response has norm bounded by (1), hence so does its restriction to \(I_i\); blocks after \(i\) do not contribute. The sum over earlier blocks is retained in full. There is no restart or deletion of old error.

Define the explicit finite constants
\[
L=2c\max_i\sum_{j\le i}\frac{w_i}{w_j}V_j,
\qquad
C=c\max_i\sum_{j\le i}\frac{w_i}{w_j^2}.                        \tag{9}
\]
Then
\[
\|2\mathcal Q(v_0,e)\|_*\le L\|e\|_*,
\qquad
\|\mathcal Q(e,f)\|_*\le C\|e\|_*\|f\|_*.
\]
Let
\[
d=\|\mathcal Hr_0\|_*.
\]
If
\[
\boxed{L<1,\qquad 4Cd<(1-L)^2,}                               \tag{10}
\]
put
\[
\rho_*=\frac{1-L-\sqrt{(1-L)^2-4Cd}}{2C},\qquad
q_*=L+2C\rho_*<1.
\]
The same invariant-ball proof now applies in \(\|\cdot\|_*\). It proves convergence of the **unpreconditioned global iteration (2)**, not merely existence of some solution constructed separately on each time interval. In particular,
\[
\boxed{\|v_{j+1}-v_j\|_*\le q_*^j d,\qquad
\|u-v_j\|_{\mathcal X}\le
K_w\,\frac{q_*^j d}{1-q_*}.}                                  \tag{11}
\]
Every \(v_j\) is a single global reference, retaining all earlier forcing through its exact causal formula. A finite partition also has no finite accumulation of endpoints.

The block defect \(d\) may be bounded using full residual-memory or signed heat-defect certificates. It is not necessary to replace it by a potentially much larger critical \(L^1\) forcing budget.

## 4. Explicit weights for any finite reference budget

Suppose \(\|v_0\|_{\mathcal X}^2=P<\infty\). Fix \(0<\alpha<1\) and \(\eta>0\) with
\[
\ell:=\frac{2c\eta}{1-\alpha}<1.
\]
Choose a finite partition such that \(V_i\le\eta\), with for example
\[
m=\max\{1,\lceil P/\eta^2\rceil\},\qquad w_i=\alpha^{i-1}.
\]
Absolute continuity of the cumulative reference integral permits such a partition; the last block may extend to infinity. Finite geometric sums in (9) give
\[
L\le\ell,\qquad C\le\overline C:=
\frac{c}{\alpha^{m-1}(1-\alpha^2)}.                            \tag{12}
\]
Consequently a fully explicit sufficient condition is
\[
\boxed{d<
\frac{(1-\ell)^2\alpha^{m-1}(1-\alpha^2)}{4c}.}                \tag{13}
\]
Because \(w_i\le1\), one may conservatively replace \(d\) by the full norm \(\|\mathcal Hr_0\|_{\mathcal X}\).

For a simple concrete choice, let \(\alpha=1/4\), \(\eta=1/\sqrt2\). Then \(\ell=2/3\), \(m=\max\{1,\lceil2P\rceil\}\), and
\[
\boxed{\|\mathcal Hr_0\|_{\mathcal X}
<\frac{5}{192c}\,4^{-(m-1)}}                                  \tag{14}
\]
guarantees convergence. These constants are conservative; the actual finite sums (9) can improve them.

Thus every specified reference with finite \(P\) has a positive certified defect neighborhood for this iteration. The allowed defect can become exponentially small as \(P\) grows. This is a local stability theorem, not evidence that an arbitrary Galerkin reference automatically has a sufficiently small defect.

## 5. The initial datum must be handled explicitly

The iteration (2) keeps the initial datum \(a\) exactly. If the desired datum is instead \(u_0=a+h\), iterating (2) alone converges, when it converges, to the solution for \(a\). It does not recover the missing initial tail.

For the target datum \(a+h\), use the map
\[
T_h(f)=S(a+h)-\mathcal Q(f,f).
\]
Relative to the same reference \(v_0\), its signed defect is
\[
g_h=T_h(v_0)-v_0=Sh-\mathcal Hr_0.
\]
The theorem holds with
\[
d_h=\|g_h\|_*
\]
in place of \(d\), with the same \(L,C\). A directly certified signed norm can retain cancellation. A simple positive majorant is
\[
d_h\le d+\|Sh\|_*\le d+\mathfrak h(h),\qquad
\mathfrak h(h)=\|Sh\|_{\mathcal X}.
\]
Hence a convenient perturbation condition is
\[
\boxed{d+\mathfrak h(h)<\frac{(1-L)^2}{4C}.}                   \tag{15}
\]
The critical bound \(\mathfrak h(h)\le A_{-1}(h)/\sqrt2\) converts this to a full-support critical perturbation radius if desired.

There are two honest implementation choices:

- Start from the auxiliary center \(v_0\) and apply \(T_h\). The first image has initial datum \(a+h\); every later residual correction preserves that new datum.
- To preserve \(a+h\) from the first iterate, start from \(\widetilde v_0=v_0+Sh\). Under the positive majorant (15), use the invariant radius obtained with \(d+\mathfrak h(h)\). It is at least \(\mathfrak h(h)\), so \(\widetilde v_0\) lies inside the same certified ball. Its subsequent steps are ordinary residual heat corrections and all have initial datum \(a+h\).

The second construction may have infinite spatial support if \(h\) does; the theorem allows this analytic seed, but it must not be described as a finite Fourier computation. A finite reference with the initial tail handled by the certificate is usually sufficient.

## 6. Explicit NS-B consequence, including every solution cutoff

Assume \(v_0\) has support \(|k|\le K\), the target \(u_0\) is smooth, and one of the target-data gates above passes with radius \(\rho_*\). Let
\[
\mathcal V=\|v_0\|_{\mathcal X}+K_w\rho_*.
\]
The fixed point satisfies \(\|u\|_{\mathcal X}\le\mathcal V\). The same bounds hold for every solution Galerkin cutoff \(M\ge K\): use the projected map \(SP_Mu_0-P_M\mathcal Q\) and the same finite reference. Projection contracts the defect and all Fourier-amplitude estimates, so every such map obeys the identical invariant-ball argument. This supplies a noncircular uniform construction.

The previously proved H¹ estimate gives
\[
\sup_t\|\nabla u_M(t)\|_2^2\le X_0e^{\mathcal V^2},
\qquad
\int_0^\infty\|\Delta u_M(t)\|_2^2dt\le X_0e^{\mathcal V^2},
\]
where \(X_0=\|\nabla u_0\|_2^2\). Combining with \(\int\|\nabla u_M\|_2^2\le E_0/2\), \(E_0=\|u_0\|_2^2\), yields
\[
\boxed{\sup_{M\ge K}\int_0^\infty\|\nabla u_M(t)\|_2^4dt
\le\frac{E_0X_0}{2}e^{\mathcal V^2}.}                         \tag{16}
\]
For the finitely many smaller cutoffs, the direct energy estimate gives \(K^2E_0^2/2\). Taking the maximum includes every cutoff. Standard strong compactness and H¹ continuation then give a global strong solution, smooth for smooth initial data. The reference and iteration estimates do not assume a global PDE solution beforehand.

## Scope and practical meaning

The basic gate (5) is the existing heat-space reference gate interpreted as an iteration theorem. The causal block theorem (10) additionally certifies convergence of the same global correction iteration around a reference with any finite \(\mathcal X\) norm, provided its complete defect is small enough. Its off-diagonal terms explicitly retain past forcing and low-frequency output.

This supplies a valid stopping/success certificate after any finite correction stage: recenter at that stage, compute or bound its complete heat defect and block reference norms, and check (10). Passing implies convergence with the explicit error estimate (11) and NS-B bound (16). Failure only means this sufficient test did not certify success.

Neither the small-defect hypothesis nor a uniform bound on the required number of corrections is proved for general initial data. Repeated correction is not automatically a convergent algorithm for arbitrary three-dimensional data. No improvement of the previously proved strict logarithmic growth threshold six is asserted here.
