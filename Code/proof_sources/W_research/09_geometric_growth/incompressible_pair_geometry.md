# Input incompressibility sharpens the symmetric Leray-pair estimate

Private written derivation. No old script, numerical integration, external source, upload, or original-file modification is used. The setting is viscosity one on the normalized periodic 2π torus. Fourier coefficients are complex vectors with the Euclidean norm; wavevectors and orthogonal projectors are real. The result concerns the same unforced, mean-zero, divergence-free NS-B problem. No novelty or arbitrary-data regularity claim is made.

## 1. Precise pair statement

Let p,q∈ℝ³ be nonzero, k=p+q≠0, and let a,b∈ℂ³ satisfy p·a=0 and q·b=0. Write P=|p|, Q=|q|, K=|k|, and let P_k be orthogonal projection onto k⊥. Define

\[
F=P_k[(k\cdot a)b+(k\cdot b)a].
\]

Then

\[
\boxed{K|F|\le(P^2+Q^2)|a||b|.}                 \tag{1}
\]

Together with the previously proved unconstrained pair bound

\[
|F|\le K|a||b|,                                  \tag{2}
\]

this improves the multiplier estimate at one derivative. Condition p·a=q·b=0 is essential to the proof of (1); (2) alone does not imply (1).

## 2. Exact complex-vector reduction

The expression is symmetric under (p,a)↔(q,b), so assume P≥Q. If p,q are collinear, both k·a and k·b vanish, and (1) is immediate. Otherwise take a real orthonormal frame (n,m,z), with n=k/K, p·m=h>0, q·m=−h, and z normal to the plane of p,q. Write

\[
p=p_n n+hm,\qquad q=q_n n-hm.
\]

The unit in-plane vectors perpendicular to p and q are

\[
e_p=(-hn+p_nm)/P,\qquad e_q=(hn+q_nm)/Q.
\]

Every permitted complex vector has the form

\[
a=a_1e_p+a_2z,\qquad b=b_1e_q+b_2z,
\]

where a₁,a₂,b₁,b₂ are complex scalars and

\[
|a|^2=|a_1|^2+|a_2|^2,\qquad
|b|^2=|b_1|^2+|b_2|^2.
\]

Let c=(p·q)/(PQ), σ=√(1−c²)>0. The exact geometric identities are

\[
Kh=PQ\sigma,\qquad p_n-q_n=(P^2-Q^2)/K=:d.
\]

Taking k·a=−Kh a₁/P and k·b=Kh b₁/Q, expanding, and deleting the n-component gives the exact formula

\[
F=\sigma\left[d a_1b_1m+(P a_2b_1-Q a_1b_2)z\right].           \tag{3}
\]

Since m,z are real orthonormal,

\[
|F|^2=\sigma^2\left[d^2|a_1b_1|^2+
                         |P a_2b_1-Q a_1b_2|^2\right].        \tag{4}
\]

Complex Cauchy–Schwarz gives

\[
|P a_2b_1-Q a_1b_2|^2
\le(P^2|a_2|^2+Q^2|a_1|^2)(|b_1|^2+|b_2|^2).
\]

Therefore

\[
\boxed{|F|^2\le\sigma^2\max\{P^2,Q^2+d^2\}|a|^2|b|^2.}       \tag{5}
\]

No restriction to real polarizations, aligned phases, or one vector component was made.

## 3. Obtuse and acute cases prove (1)

Put r=Q/P∈(0,1].

**If c≤0**, the already proved bound (2) immediately gives

\[
K|F|\le K^2|a||b|
=(P^2+Q^2+2PQc)|a||b|
\le(P^2+Q^2)|a||b|.                             \tag{6}
\]

**If c≥0**, then K²≥P²+Q²≥P²−Q², so d²≤P²−Q² and the maximum in (5) is P². It remains to bound KPσ. The exact nonnegative identity

\[
\begin{aligned}
&(1+r^2)^2-(1+r^2+2rc)(1-c^2)\\
&\hspace{1em}=(r^2+c^2)^2+(r-c)^2(1-c^2)\ge0
\end{aligned}                                                     \tag{7}
\]

shows that K²P²σ²≤(P²+Q²)². Multiplying (5) by K² proves (1) also in the acute case. These two cases include c=0 consistently.

The collinear cases and k=0 cause no hidden division. The collinear cases were treated first; the zero output mode of the NS nonlinearity is identically zero and is omitted from every multiplier sum.

The numerical coefficient 1 in (1) cannot be reduced as a uniform pair inequality. For p=(n,0,0), q=(0,1,0), a=(0,0,1), b=(1,0,0), one has F=na, so

\[
\frac{K|F|}{(P^2+Q^2)|a||b|}
=\frac{n}{\sqrt{n^2+1}}\longrightarrow1.
\]

This proves sharpness of this pair coefficient; it is not a separate claim of sharpness for every subsequent convolution or energy estimate.

## 4. Consequence for the quadratic NS multiplier

Let v be any smooth or finite Fourier divergence-free field, and put U=A₀(v)=Σ_{k≠0}|v̂_k|. Symmetrizing the ordered Fourier convolution gives

\[
\widehat{B(v,v)}_k
=\frac i2\sum_{p+q=k}P_k[(k\cdot\widehat v_p)\widehat v_q
                         +(k\cdot\widehat v_q)\widehat v_p].
\]

By (1),

\[
|k|\,|\widehat{B(v,v)}_k|
\le\frac12\sum_{p+q=k}(|p|^2+|q|^2)|\widehat v_p||\widehat v_q|
=\sum_{p+q=k}|p|^2|\widehat v_p||\widehat v_q|.
\]

Discrete Young ℓ²*ℓ¹→ℓ² now proves

\[
\boxed{\|\Lambda B(v,v)\|_2
       \le A_0(v)\|\Lambda^2v\|_2,\qquad\Lambda=(-\Delta)^{1/2}.} \tag{8}
\]

The former safe coefficient 2 is replaced by 1. The pair symmetrization factor 1/2 is retained, and the two weighted convolution terms combine to exactly one copy. No extra factor 2 remains.

More generally, for 0<s≤2, (1)–(2) imply

\[
\begin{aligned}
K^{s-1}|F|
&\le\min\{K^s,(P^2+Q^2)K^{s-2}\}|a||b|\\
&\le(P^2+Q^2)^{s/2}|a||b|\\
&\le(P^s+Q^s)|a||b|.
\end{aligned}                                                     \tag{9}
\]

For the middle inequality, split at K=√(P²+Q²): use the first term below that point and the second term above it. The last inequality is concavity of x↦x^(s/2). Consequently

\[
\boxed{\|\Lambda^{s-1}B(v,v)\|_2
       \le A_0(v)\|\Lambda^s v\|_2\quad(0<s\le2).}               \tag{10}
\]

At s≥2, a safe extension is C_s=2^(s−2), using (1), K≤P+Q, and finite-dimensional norm comparison:

\[
(P+Q)^{s-2}(P^2+Q^2)\le2^{s-2}(P^s+Q^s).
\]

No optimal constant for s>2 is asserted or needed here. This extension has no better cutoff threshold than s=2.

## 5. Improved H² weighted residual gate

Take the exact spherical Galerkin reference v=v_N with v(0)=P_Nu₀. Define

\[
Y=\|\Delta v\|_2^2,\quad Z=\|\Lambda^3v\|_2^2,\quad
U=A_0(v),\quad a(t)=\int_0^tU^2,\quad\mathcal A_N=a(\infty),
\quad Y_0=\|\Delta u_0\|_2^2.
\]

For each finite N these reference quantities are well defined, and the exact Galerkin energy equality gives finite 𝒜_N. Testing at H² and using (8) gives

\[
Y'+2Z\le2U\sqrt{YZ},
\]

so for any 0<η<2,

\[
Y'+(2-\eta)Z\le\eta^{-1}U^2Y.                    \tag{11}
\]

For the full truncation residual r_N=(I−P_N)B(v,v), let R_N=A₋₁(r_N). Every output N<|k|≤2N is included. At least one interacting retained input lies in N/2<|p|≤N, hence

\[
R_N\le U\sqrt{S_6^{\rm band}(N)}\sqrt Z,
\quad S_6^{\rm band}(N)=\sum_{N/2<|k|\le N}|k|^{-6}.             \tag{12}
\]

The full initial tail satisfies

\[
b_0=A_{-1}((I-P_N)u_0)
\le\sqrt{Y_0S_6^{\rm tail}(N)},\quad
S_6^{\rm tail}(N)=\sum_{|k|>N}|k|^{-6}.                         \tag{13}
\]

For every β>1/4 choose α=√β and η=1/√β. Then η<2, α>1/2, and 2β−α>0. Multiplying (11) by exp(−αa), integrating, and applying time Cauchy–Schwarz to (12) proves

\[
\boxed{
F_{\beta,N}:=b_0+\int_0^\infty e^{-\beta a(t)}R_N(t)\,dt
\le\sqrt{Y_0}\left[
\sqrt{S_6^{\rm tail}(N)}+
\frac{\sqrt{S_6^{\rm band}(N)}}{2\sqrt\beta-1}\right].}         \tag{14}
\]

Specifically the denominator before optimization is

\[
\sqrt{(2-1/\alpha)(2\beta-\alpha)},
\]

whose maximum for α∈(1/2,2β) is 2√β−1 at α=√β. Thus every sign, admissibility condition, and constant in (14) is explicit.

The existing critical error inequality is

\[
b'+(1-\varepsilon-b/2)h
\le\frac{U^2}{4\varepsilon}b+R_N,
\quad b=A_{-1}(e),\quad h=A_1(e).
\]

With ε=1/(4β)∈(0,1), its first-contact argument closes whenever

\[
\boxed{
e^{\beta\mathcal A_N}\sqrt{Y_0}
\left[\sqrt{S_6^{\rm tail}(N)}+
\frac{\sqrt{S_6^{\rm band}(N)}}{2\sqrt\beta-1}\right]
<2-\frac1{2\beta},\qquad\beta>\frac14.}                         \tag{15}
\]

A smaller independently known b₀ can replace the majorant in (13). The error inequality applies to all higher Galerkin approximations against this same reference with contracted initial tail and residual. Thus the continuation conclusion does not assume a global PDE solution in advance. For smooth data, a passing gate proves global smooth regularity; for H² data it gives the corresponding global strong solution and positive-time regularity.

For every real N≥1 the previously established elementary constants give

\[
S_6^{\rm tail}(N)\le\frac{1000}{7N^3},\qquad
S_6^{\rm band}(N)\le\frac{1728}{N^3}.
\]

Hence one entirely explicit sufficient inequality is

\[
e^{\beta\mathcal A_N}\sqrt{Y_0}N^{-3/2}
\left[10\sqrt{\frac{10}{7}}+
\frac{24\sqrt3}{2\sqrt\beta-1}\right]
<2-\frac1{2\beta}.                                               \tag{16}
\]

## 6. Cutoff threshold, unweighted defect, and exact scope

For a fixed H² datum, (15) proves the stronger sufficient condition

\[
\boxed{\liminf_{N\to\infty}
       \frac{\mathcal A_N}{\log N}<6.}                           \tag{17}
\]

Indeed choose a positive ℓ strictly between the liminf and 6, then choose fixed β with 1/4<β<3/(2ℓ). Along a subsequence 𝒜_N≤ℓ log N, and the left side of (16) tends to zero as N^(βℓ−3/2). The right side is a fixed positive number. A limsup below 6 gives passage at all sufficiently large cutoffs. The zero datum is trivial. Equality at 6 is not asserted: β↓1/4 makes both the denominator in the residual budget singular and the allowable error barrier tend to zero.

Within the safe constants proved here, (10) gives threshold 4(s−1/2) for 1/2<s≤2, maximized at 6. The extension C_s=2^(s−2) gives 4^(3−s)(s−1/2) for s≥2, decreasing there. Thus 6 is the maximum of this proved family, not a claim of a sharp PDE threshold.

The same improvement also strengthens the honest unweighted heat-space defect. With

\[
D_N=\left\|\int_0^t e^{(t-\tau)\Delta}r_N(\tau)\,d\tau\right\|_{L_t^2A_0},
\]

modewise heat convolution and (11) give, for every α>1/2,

\[
D_N^2\le\frac{Y_0S_6^{\rm band}(N)}2
                \frac{e^{\alpha\mathcal A_N}-1}{2\alpha-1}.     \tag{18}
\]

For 𝒜_N>0, choosing α=1/2+1/(2𝒜_N) yields

\[
D_N^2\le\frac{Y_0S_6^{\rm band}(N)\mathcal A_N}{2}
        \left[e^{\mathcal A_N/2+1/2}-1\right].                  \tag{19}
\]

If 𝒜_N=0, the reference and residual vanish. Thus the unweighted defect has a bound of order N^(−3/2)√(Y₀𝒜_N)exp(𝒜_N/4). This retains amplification; it is not an unconditional cutoff-decaying bound after its 𝒜_N dependence is dropped.

A passing gate also supplies an explicit bound for the original squared-enstrophy quantity. Let M be the left side of (15), including a chosen certified upper bound for its forcing budget, and define

\[
\varepsilon=1/(4\beta),\qquad
H_*={M\over1-\varepsilon-M/2},\qquad
\Gamma_N=\left(\sum_{0<|k|\le N}|k|^2\right)^{1/2}.
\]

The gate guarantees a positive denominator. Its critical error proof gives b≤M and ∫A₁(e)≤H_*. Exact energy decay of the reference gives

\[
\int_0^\infty A_1(v_N)\,dt\le\Gamma_N\sqrt{E_0}.
\]

Consequently, for every Galerkin cutoff m≥N, the H¹ inequality

\[
X_m'+2Y_m\le2A_1(u_m)X_m
\]

gives sup X_m≤X₀ exp[2(Γ_N√E₀+H_*)]. For m<N the direct reference estimate ∫A₁(u_m)≤Γ_m√E₀≤Γ_N√E₀ gives the same upper bound without the error term. Finally ∫X_m≤E₀/2 for every m, so

\[
\boxed{\sup_{m\ge1}\int_0^\infty\|\nabla u_m(t)\|_2^4\,dt
\le{E_0X_0\over2}
\exp\!\left[2\Gamma_N\sqrt{E_0}+2H_*\right],
\quad X_0=\|\nabla u_0\|_2^2.}                   \tag{20}
\]

This bound may be conservative, but it is finite, explicit after a passing certificate, and independent of the solution cutoff m. Its dependence on the selected passing reference N is allowed. It preserves the original NS-B objective rather than substituting only a regularity statement.

The unconditional energy estimate is still only

\[
\mathcal A_N\le\tfrac12E_0S_2(N)=O(E_0N).
\]

No estimate proving (17) for every initial datum has been obtained. Common-time eventual entry of finite references also does not bound their preceding accumulated 𝒜_N uniformly in N. The completed progress is the rigorous incompressibility improvement (1), its nonlinear norm consequence (8), and the stronger sufficient general-data bridge (15)–(17). Universal NS-B remains open.
