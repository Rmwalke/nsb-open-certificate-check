# One exact heat-defect correction of a Galerkin reference

Private written NS-B continuation. No old research script, old import, numerical time integration, external search, upload, or original-file modification is used. One fresh finite rational calculation checked the supplementary example below. The setting is viscosity one on the normalized 2π torus, with real mean-zero divergence-free fields and Euclidean complex Fourier coefficient norms.

The correction can improve a concrete certificate. Its defect is not automatically small, and its new residual need not remain at high frequencies.

## 1. Exact equation, support, and initial error

Let \(v=v_N\) be the exact spherical Galerkin solution with initial datum \(P_Nu_0\), and use the unsymmetrized notation \(B(a,b)=\mathbb P((a\cdot\nabla)b)\). With \(A=-\Delta\), define
\[
r=v'+Av+B(v,v)=(I-P_N)B(v,v),\qquad
z'+Az=r,\quad z(0)=0,\qquad w=v-z.
\]
Direct substitution, without dropping any mode, proves
\[
\boxed{r_w:=w'+Aw+B(w,w)=-B(v,z)-B(z,v)+B(z,z).}       \tag{1}
\]

| Quantity | Fourier support |
|---|---|
| \(v\) | \(0<|k|\le N\) |
| \(r,z\) | \(N<|k|\le2N\) |
| \(w\) | \(0<|k|\le2N\) |
| \(r_w\) | \(0<|k|\le4N\), potentially including low modes |

Zero output modes vanish by divergence form. These are upper support bounds; actual supports may be smaller. Since \(z(0)=0\),
\[
\boxed{u_0-w(0)=(I-P_N)u_0.}                         \tag{2}
\]
The initial truncation error is unchanged. Restarting a fresh correction at a later time changes the reference and must include the old correction in its reference jump; it cannot silently discard preceding error.

## 2. Reference amplification versus forcing reduction

Write
\[
U=A_0(v),\quad \zeta=A_0(z),\quad
P=\int_0^\infty U^2,\quad D=\|z\|_{L_t^2A_0},\quad
U_w=A_0(w),\quad P_w=\int_0^\infty U_w^2.
\]
The supports of \(v,z\) are disjoint, so these identities are exact:
\[
\boxed{U_w=U+\zeta,\qquad
P_w=P+2\int U\zeta+D^2,\qquad
P\le P_w\le(\sqrt P+D)^2.}                         \tag{3}
\]
The same disjoint-support identity holds for \(A_{-1}\) and \(A_1\). Thus no Fourier-amplitude cancellation between \(v\) and \(-z\) reduces the reference growth budget.

The symmetric Leray-pair estimate gives
\[
\boxed{R_w:=A_{-1}(r_w)\le U\zeta+\tfrac12\zeta^2
       =\zeta(U+\zeta/2).}                         \tag{4}
\]
In particular,
\[
\boxed{\int R_w\le K_w:=\int U\zeta+\tfrac12D^2
=\tfrac12(P_w-P)\le\sqrt P\,D+\tfrac12D^2.}        \tag{5}
\]
This is an estimate about a specified reference, not an assumption that \(D\) is small. Under small-amplitude scaling the original residual starts at quadratic order and the corrected residual at cubic order. At general amplitude, (3)–(5) record the competing costs honestly.

## 3. Complete causal budgets for the correction

Put \(R=A_{-1}(r)\), \(\lambda=N^2\), and
\[
S_b=\sum_{N/2<|p|\le N}|p|^{-6},\qquad
Z=\|\Lambda^3v\|_2^2,\quad \Lambda=(-\Delta)^{1/2}.
\]
Every original residual output is retained, and its input localization is
\[
R\le U\sqrt{S_b}\sqrt Z.
\]
The band is an interior retained input band, not the empty set of reference inputs above \(N\). The causal residual-memory theorem gives
\[
A_{-1}(z(T))\le C(T):=\int_0^T e^{-N^2(T-s)}R(s)\,ds,
\qquad \mathfrak h(z(T))\le C(T)/\sqrt2,           \tag{6}
\]
and
\[
\boxed{D^2\le\frac12\iint e^{-N^2|s-r|}R(s)R(r)\,ds\,dr
\le\min\{N^{-2}\|R\|_{L^2}^2,\tfrac12\|R\|_{L^1}^2\}.}        \tag{7}
\]
Here \(\mathfrak h(f)=\|e^{t\Delta}f\|_{L_t^2A_0}\). The double integral follows by squaring the positive Fourier heat majorant and integrating its individual kernel overlaps; no invalid Fourier/time norm exchange is used.

For \(Y_0=\|\Delta u_0\|_2^2\), the coefficient-one H² estimate and completed-square covariance also imply
\[
\boxed{D\le\frac3{\sqrt2}\sqrt{Y_0S_b}\,e^{P/4}.}              \tag{8}
\]
The sharper finite-window covariance in causal_residual_memory.md may replace this coarse bound. Equations (6)–(8) make \(z\) a completely specified correction; they do not make it cutoff-small when \(P\) is uncontrolled.

## 4. A causal weighted forcing gain

Define
\[
a(t)=\int_0^tU^2,\qquad a_w(t)=\int_0^t(U+\zeta)^2.
\]
For each original residual frequency \(k\), \(N<|k|\le2N\), and \(s<T\le\infty\), put
\[
G_{k,T}(s)=\int_s^T e^{-[a_w(t)-a_w(s)]/4}
(U(t)+\zeta(t)/2)|k|e^{-|k|^2(t-s)}\,dt,
\]
and let \(G_T(s)\) be its supremum over those frequencies. Substituting the positive heat majorant of \(\zeta\) into (4), and writing \(f_k=|\widehat r_k|/|k|\), gives by Tonelli
\[
\boxed{\int_0^T e^{-a_w(t)/4}R_w(t)\,dt
\le\int_0^T e^{-a_w(s)/4}G_T(s)R(s)\,ds.}          \tag{9}
\]
This keeps every forcing time and its subsequent interaction with the reference. Because \(U+\zeta/2\le U+\zeta\) and \(a_w'=(U+\zeta)^2\), time Cauchy–Schwarz proves
\[
\begin{aligned}
G_{k,T}(s)^2
&\le\left[\int_s^T e^{-[a_w(t)-a_w(s)]/2}(U+\zeta/2)^2dt\right]
\left[\int_s^T |k|^2e^{-2|k|^2(t-s)}dt\right]\\
&\le(1-e^{-[a_w(T)-a_w(s)]/2})(1-e^{-2|k|^2(T-s)})
\le1-e^{-[a_w(T)-a_w(s)]/2}.
\end{aligned}                                                   \tag{10}
\]
Thus \(G_\infty(s)<1\) for finite \(P_w\). Independently,
\[
G_{k,T}(s)\le
\frac{\sup_{s\le t\le T}(U(t)+\zeta(t)/2)}{|k|}
\le\frac{\sup_{s\le t\le T}(U(t)+\zeta(t)/2)}N.       \tag{11}
\]
Actual causal profiles can improve these envelopes. The high-frequency gap here belongs to \(z\)'s response, not to \(r_w\).

Since \(a_w\ge a\), a rigorous comparison of the weighted forcing alone is
\[
\boxed{\int_0^\infty e^{-a_w/4}R_w
\le\int_0^\infty e^{-a_w/4}G_\infty R
\le\int_0^\infty e^{-a/4}R.}                       \tag{12}
\]
The complete certificate need not improve monotonically: its prefactor changes from \(e^{P/4}\) to \(e^{P_w/4}\), and the unchanged initial error is also multiplied by that larger factor.

## 5. Critical gate and the original NS-B bound

Let \(b_0=A_{-1}((I-P_N)u_0)\) and \(\Psi(b)=be^{-b/2}\). The endpoint critical comparison proof applies to any smooth finite reference with its full residual; exact Galerkin dynamics of that reference are unnecessary. It therefore applies to \(w\) under
\[
\boxed{b_0<2,\qquad
M_w:=e^{P_w/4}\left[\Psi(b_0)+\int_0^\infty e^{-a_w(t)/4}R_w(t)dt\right]
<2/e.}                                                        \tag{13}
\]
The residual may be enclosed directly from (1), by (4), or by the causal gain (9). The initial condition \(b_0<2\) is retained because \(\Psi\) is not increasing beyond two.

Let \(B=\Psi^{-1}(M_w)\) on its increasing branch \([0,2)\), \(\delta=1-B/2\), and
\[
H_*=\frac{2b_0}{\delta}+\frac{BP_w}{\delta^2}+\frac{2K_w}{\delta}.
\]
The critical error proof gives \(b(t)\le B\), \(\int A_1(u-w)\le H_*\), and
\[
\|u\|_{L_t^2A_0}\le\sqrt{P_w}+\sqrt{BH_*}.
\]
These inequalities first hold uniformly for every solution cutoff \(m\ge2N\), by projecting the same residual (1), and then pass to the PDE. With \(E_0=\|u_0\|_2^2\), \(X_0=\|\nabla u_0\|_2^2\), they yield
\[
\boxed{\sup_{m\ge2N}\int_0^\infty\|\nabla u_m(t)\|_2^4dt
\le\frac{E_0X_0}{2}\exp[(\sqrt{P_w}+\sqrt{BH_*})^2].}          \tag{14}
\]
For \(m<2N\), the direct Fourier/energy estimate gives \(\int\|\nabla u_m\|_2^4\le2N^2E_0^2\). Taking the maximum includes every cutoff. Thus a passing gate gives an actual finite NS-B bound and smooth-data global continuation.

## 6. Supplementary exact example: the same gate fails before and passes after correction

This example uses coarse pair estimates, rather than exact cubic and quartic residual coefficients. It demonstrates improvement of this particular endpoint criterion; it does not claim that every earlier or higher-cutoff certificate fails.

Let
\[
W=(\sin z,\sin x,\sin y),\quad
T=(\sin y\cos z,\sin z\cos x,\sin x\cos y),\quad a=37/50.
\]
Direct differentiation gives \(B(W,W)=T\). The latter is already divergence-free, has all frequencies of squared length two, and \(A_0(W)=A_0(T)=3\). Since \(P_1T=0\), the exact cutoff-one reference and heat correction are
\[
v(t)=ae^{-t}W,\quad r(t)=a^2e^{-2t}T,\quad
z(t)=a^2t e^{-2t}T.                                \tag{15}
\]
The initial datum \(u_0=aW\) has \(b_0=0\). All formulas are closed analytic expressions; no ODE integration is used. Exact integration gives
\[
P=\frac92a^2,\qquad P_w=P+\Delta P,\qquad
\Delta P=2a^3+\frac9{32}a^4.                        \tag{16}
\]

### Original gate fails

Here \(R=(3/\sqrt2)a^2e^{-2t}\), \(a(t)=P(1-e^{-2t})\), and the exact original gate is
\[
M_v=\frac4{3\sqrt2}(e^{9a^2/8}-1).
\]
Put \(x=9a^2/8=12321/20000\). Exact arithmetic gives
\[
x+x^2/2=\frac{644647041}{800000000}>\frac{45}{56}.
\]
Since \(\sqrt2<10/7\) and \(e>8/3\),
\[
\frac{3\sqrt2}{2e}<\frac{45}{56}<x+x^2/2<e^x-1.
\]
Thus \(M_v>2/e\), so this specific uncorrected endpoint gate fails.

### Corrected gate passes

Equation (4) gives
\[
R_w(t)\le9a^3t e^{-3t}+\frac92a^4t^2e^{-4t}.
\]
Using \(a_w(t)\ge a(t)\), expanding the positive exponential, and applying Tonelli,
\[
M_w\le e^{\Delta P/4}\sum_{n=0}^\infty\frac{x^n}{n!}
\left[\frac{9a^3}{(3+2n)^2}+\frac{9a^4}{(4+2n)^3}\right].     \tag{17}
\]
The integrations used are exactly \(\int t e^{-ct}=c^{-2}\) and \(\int t^2e^{-ct}=2c^{-3}\).

Rational majorants suffice:
\[
x<5/8,\quad a^3<203/500,\quad a^4<3/10,\quad
\Delta P=\frac{178957049}{200000000}<1.
\]
Define
\[
U_n=\frac{(5/8)^n}{n!}
\left[\frac{1827/500}{(3+2n)^2}+\frac{27/10}{(4+2n)^3}\right].
\]
For \(n\ge2\), monotonicity of the denominators gives \(U_{n+1}/U_n\le5/24\). Hence
\[
\sum U_n\le U_0+U_1+\frac{24}{19}U_2
=\frac{772274787}{1361920000}<\frac47.              \tag{18}
\]
Similarly,
\[
e^{5/4}\le\sum_{j=0}^4\frac{(5/4)^j}{j!}
+\frac{24}{19}\frac{(5/4)^5}{5!}
=\frac{135827}{38912}<\frac72,                     \tag{19}
\]
because after the fifth exponential term the successive ratio is at most \(5/24\). Therefore
\[
\boxed{eM_w<e^{5/4}\sum U_n<\frac72\frac47=2.}       \tag{20}
\]
The product of the two explicit rational majorants is
\[
\frac{104895767493849}{52995031040000}<2.
\]
Thus the corrected reference passes the same endpoint gate that failed for the original exact reference. This is a real finite certificate improvement for \((37/50)W\), not merely a formal rewriting of the mild equation.

### Its new residual actually contains low modes

Direct differentiation gives
\[
P_1[B(W,T)+B(T,W)]=-\tfrac12W.
\]
For example, the first component before projection is \(2\sin x\cos y\cos z-\sin^2y\sin z\), whose unit-frequency part is \(-\sin z/2\); the other components cycle. These unit-frequency vector modes are already divergence-free. Meanwhile \(B(T,T)\) has even coordinate-sum parity and no unit-axis frequency. Consequently
\[
\boxed{P_1r_w(t)=\tfrac12a^3t e^{-3t}W\ne0.}        \tag{21}
\]
The original residual's high-output gap therefore cannot be assigned to the corrected residual.

## 7. Worst-case growth and the existing threshold six

For a fixed H² datum, suppose along an unbounded cutoff subsequence \(P_N\le\ell\log N\) for some fixed \(\ell<6\). The original endpoint budget is at most a constant times \(e^{P_N/4}N^{-3/2}(1+P_N)\), which tends to zero. Equation (8) gives
\[
D_N=O(N^{-3/2}e^{P_N/4})\to0,\qquad
P_{w,N}-P_N\le2\sqrt{P_N}D_N+D_N^2\to0.
\]
From (12), with the unchanged initial term,
\[
M_{w,N}\le e^{(P_{w,N}-P_N)/4}M_{v,N}\longrightarrow0.
\]
Thus the correction preserves the sufficient strict slope-six bridge, and the initial \(b_0<2\) is eventually automatic from the H² tail bound.

These estimates do not improve that universal sufficient slope: the guaranteed gain \(G\le1\) supplies no extra cutoff-decaying power. The stronger bound (11) requires a reference amplitude envelope not uniformly supplied by energy. Using only \(P_N=O(E_0N)\) in (8) is too weak even to prove that \(D_N\) is small. No automatic success, universal convergence of repeated corrections, or nonaccumulation theorem follows.

The completed result is an exact finite correction with full support and initial-error accounting, causal and weighted forcing estimates, an NS-B gate, and a concrete strict certificate improvement. Universal NS-B remains open.
