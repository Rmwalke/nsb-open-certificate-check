# An analytic residual-corrected certificate for W and a critical neighborhood

Private written result, 2026-09-30 (Asia/Kolkata). The derivation uses
fresh exact Fourier algebra and finite rational series bounds. It uses
no numerical ODE/PDE evolution, no previous research-code execution,
no external source, and no previous W numerical certificate.

The endpoint gate passes for the corrected reference at amplitude one.
In fact it gives a uniform analytic certificate for

\[
 \boxed{u_0=aW+h,\quad |a|\le1,\quad
 A_{-1}(h)\le1/500,\qquad W=(\sin z,\sin x,\sin y),}           \tag{1}
\]

where \(h\) is any smooth real, mean-zero, divergence-free field.
There is no Fourier-support or symmetry restriction on \(h\).
The constants below are deliberately conservative. This is an analytic
certificate for a specified class, not an arbitrary-data result or a
novelty claim.

## 1. Full defect correction and exact residual

Use viscosity one and normalized measure on the \(2\pi\)-torus. Let
\(B(u,v)=\mathbb P((u\cdot\nabla)v)\), \(A=-\Delta\), and
\(A_s(u)=\sum_{k\ne0}|k|^s|\widehat u_k|\).

The first heat reference is \(v(t)=a e^{-t}W\). Put

\[
 F=B(W,W)=(\sin y\cos z,\sin z\cos x,\sin x\cos y).
\]

This vector is already divergence-free, \(AF=2F\), and
\(A_0(F)=3\). The full residual of \(v\) is
\(r_v=a^2e^{-2t}F\). Therefore the exact linear residual response is

\[
 z(t)=a^2t e^{-2t}F,\qquad z'+Az=r_v,\qquad z(0)=0.
\]

Define the corrected finite Fourier reference

\[
 \boxed{w(t)=a e^{-t}W-a^2t e^{-2t}F.}              \tag{2}
\]

Its initial value is \(aW\), so the complete initial error for (1) is
\(h\). Its full residual, with no omitted output modes, is

\[
 \boxed{r_w=-a^3t e^{-3t}G+a^4t^2e^{-4t}H,\quad
 G=B(W,F)+B(F,W),\quad H=B(F,F).}                   \tag{3}
\]

The linear residual cancels exactly because of the sign in (2).

## 2. Exact Fourier constants and geometric cancellation

The triple-frequency contribution to the raw mixed term is

\[
 2(\sin x\cos y\cos z,\sin y\cos z\cos x,
                         \sin z\cos x\cos y)
       =\nabla[-2\cos x\cos y\cos z].
\]

Leray projection removes it. The remaining mixed term is

\[
 G=-\tfrac12W+	frac12(
       \sin z\cos2y,\sin x\cos2z,\sin y\cos2x).
\]

It is divergence-free. Its six fundamental coefficients have magnitude
\(1/4\), and its twelve coefficients at squared frequency five have
magnitude \(1/8\). Thus

\[
 A_0(G)=3,\qquad g:=A_{-1}(G)
             =\frac32+\frac3{2\sqrt5}.
\]

The projected quadratic term \(H\) has exactly 24 nonzero modes, all
at squared frequency six. Each coefficient has squared magnitude
\(1/192\), or magnitude \(\sqrt3/24\). For example, at
\(k=(\sigma_x,\sigma_y,2\sigma_z)\), with every \(\sigma_i=\pm1\),
the projected coefficient is

\[
 -\frac{i\sigma_z}{24}
   (1,\sigma_x\sigma_y,-\sigma_x\sigma_z).
\]

The other two frequency families are obtained by cyclic permutation.
Consequently

\[
 A_0(H)=\sqrt3,\qquad q:=A_{-1}(H)=\frac1{\sqrt2}.
\]

The supports of \(G\) and \(H\) are disjoint. The exact full residual
critical norm is therefore

\[
 \boxed{R_a(t)=g|a|^3t e^{-3t}+q|a|^4t^2e^{-4t}.}              \tag{4}
\]

The fundamental and diagonal supports of \(W,F\) are also disjoint,
so the exact reference amplitude and its integral are

\[
 U_a(t)=A_0(w(t))=3|a|e^{-t}+3a^2t e^{-2t},
\]
\[
 \boxed{P_a=\int_0^\infty U_a^2dt
      =\frac92a^2+2|a|^3+\frac9{32}a^4.}           \tag{5}
\]

All these coefficient identities are independently checked by fresh
exact rational Fourier calculations. The calculation includes Leray
projection before taking coefficient norms.

## 3. Analytic endpoint residual bound at amplitude one

Let \(J_a(t)=\int_0^t U_a^2\). The endpoint critical gate uses

\[
 M_a=e^{P_a/4}\left[\Psi(b_0)+
        \int_0^\infty e^{-J_a(t)/4}R_a(t)dt\right],
 \quad \Psi(b)=be^{-b/2},\quad b_0=A_{-1}(h).
\]

It passes when \(b_0<2\) and \(M_a<2/e\). The proof then keeps the
complete error below two and controls its integrated dissipation.

First take \(a=1\). Formula (5) gives
\(P_1=217/32\). Since \(U_1(t)\ge3e^{-t}\),

\[
 J_1(t)\ge\frac92(1-e^{-2t}),
\]

and hence, with \(c=9/8\) and \(\Delta=73/128\),

\[
 e^{(P_1-J_1(t))/4}
       \le e^\Delta\exp(c e^{-2t}).
\]

Expand the last exponential into its positive series. Tonelli is
applicable, and the elementary moments
\(\int_0^\infty t e^{-\lambda t}dt=\lambda^{-2}\) and
\(\int_0^\infty t^2 e^{-\lambda t}dt=2\lambda^{-3}\) give

\[
 M_{\rm res}:=\int_0^\infty e^{(P_1-J_1(t))/4}R_1(t)dt
 \le e^\Delta\sum_{n=0}^\infty\frac{c^n}{n!}
   \left[\frac g{(3+2n)^2}+\frac{2q}{(4+2n)^3}\right].        \tag{6}
\]

This calculation integrates known exponential-polynomial functions;
it is not a numerical time evolution or an enclosure inferred from
sampled times.

Use the simple rational bounds \(g<87/40\), \(q<71/100\).
For

\[
 t_n=\frac{(9/8)^n}{n!}
       \left[\frac{87/40}{(3+2n)^2}
                   +\frac{142/100}{(4+2n)^3}\right],
\]

the successive ratio is at most \((9/8)/(n+1)\), because both
denominator factors increase with \(n\). In particular,

\[
 \sum_{n=0}^\infty t_n
 \le\sum_{n=0}^{8}t_n+
           \frac{t_9}{1-(9/8)/10}
 <\frac{51}{125}.                                 \tag{7}
\]

The final inequality is a finite exact rational comparison. A positive
Taylor series with its geometric tail also gives

\[
 e^{73/128}<177/100,\qquad
 e^{217/128}<11/2,\qquad e<68/25.
\]

For clarity, one common reproducible upper bound used for each positive
\(x<14\) is

\[
 e^x\le\sum_{n=0}^{12}\frac{x^n}{n!}
       +\frac{x^{13}/13!}{1-x/14}.
\]

Thus (6)--(7) prove

\[
 \boxed{M_{\rm res}<\frac{9027}{12500}.}            \tag{8}
\]

## 4. Uniform amplitude segment and full critical perturbations

For \(|a|\le1\), \(U_a(t)\le U_1(t)\) pointwise. Therefore the
remaining budget
\(P_a-J_a(t)=\int_t^\infty U_a^2\) is at most its value for amplitude
one. Also \(R_a(t)\le R_1(t)\). Hence the actual positive integrand
defining \(M_{\rm res}\) is monotone in \(|a|\), and (8) applies
uniformly on the entire segment.

For \(b_0\le1/500\), \(\Psi(b_0)\le b_0\), while
\(e^{P_a/4}\le e^{217/128}<11/2\). Consequently

\[
 \boxed{M_a<\frac{9027}{12500}+\frac{11}{1000}
       =\frac{18329}{25000}
       <\frac{25}{34}<\frac2e.}                   \tag{9}
\]

The separate initial condition \(b_0<2\) is automatic. This proves
the sufficient gate for every datum in (1), including perturbations at
arbitrarily high Fourier frequencies. The real and divergence-free
conditions are preserved; no symmetry of the perturbation is imposed.

This conclusion does not require a small full \(H^1\) norm of \(h\).
For example the smooth field \(h_K=(K/500)\sin(Kx_1)e_2\) satisfies
\(A_{-1}(h_K)=1/500\), while its \(H^1\) norm grows with \(K\).
The energy and enstrophy constants in the final bound below are those
of the actual full datum and may be large.

## 5. Explicit error and original NS-B constants

Choose the uniform error barrier

\[
 B_*=15/8,\qquad \delta=1-B_*/2=1/16.
\]

The same finite Taylor-tail estimate proves
\(e^{15/16}<511/200\), so

\[
 \Psi(B_*)>375/511>18329/25000.
\]

The first-contact proof of the endpoint gate therefore gives
\(A_{-1}(u-w)\le B_*\) uniformly for every higher Galerkin cutoff.
Furthermore, from (4),

\[
 \int_0^\infty R_a(t)dt
 =\frac{g|a|^3}{9}+\frac{q|a|^4}{32}
 <\frac{2533}{9600}<\frac{53}{200}.
\]

Let \(h_e=A_1(u-w)\). Return to the unsplit critical error estimate,
use the uniform barrier, and apply Young's inequality. It gives

\[
 \int_0^\infty h_e dt
 \le\frac{B_*P_a}{\delta^2}
       +\frac2\delta\left(b_0+\int_0^\infty R_a dt\right)
 <3255+32\frac{267}{1000}<3264.                    \tag{10}
\]

Thus \(\int A_0(u-w)^2\le B_*\int h_e<6120\). Minkowski in time
implies

\[
 \int_0^\infty A_0(u)^2dt
 \le(\sqrt{P_a}+\sqrt{6120})^2
 <\left(\frac{21}{8}+\frac{313}{4}\right)^2
 =\frac{418609}{64}<6541.                          \tag{11}
\]

Every reference in (2) has support inside \(|k|\le\sqrt2\), so this
comparison applies uniformly to every spherical Galerkin cutoff
\(K\ge\sqrt2\). The smaller cutoffs contain only unit frequencies
(or none). For them

\[
 A_0(P_Ku_0)=A_{-1}(P_Ku_0)
 \le3|a|+1/500<4,
 \quad \|P_Ku_0\|_2^2<16,
 \quad S_2(K)\le6.
\]

Their energy identity therefore gives
\(\int A_0(u_K)^2dt\le S_2(K)\|P_Ku_0\|_2^2/2<48\).
This argument bounds only the low-frequency projection; it does not
assert that the full datum has bounded energy throughout the critical
perturbation class.

Set \(E_0=\|u_0\|_2^2\), \(X_0=\|\nabla u_0\|_2^2\). The
direct enstrophy estimate
\(X_K'+Y_K\le A_0(u_K)^2X_K\), together with energy equality,
now gives the explicit original target

\[
 \boxed{\sup_K\int_0^\infty\|\nabla u_K(t)\|_2^4dt
              \le\frac{E_0X_0}{2}e^{6541}.}         \tag{12}
\]

The exponential constant is very large and not proposed as sharp.
Its purpose is an entirely analytic, cutoff-independent certificate.
Uniform Galerkin bounds give the global strong limit, and smooth data
give the global smooth solution. For \(u_0=W\),
\(E_0=X_0=3/2\), so (12) specializes to \((9/8)e^{6541}\).
The older validated numerical certificate has much better constants;
the present proof is independent of that computation.

## 6. What the correction improves, and what remains unresolved

For the uncorrected heat reference at \(u_0=W\), the exact endpoint
majorant is

\[
 M_{\rm heat}=\frac4{3\sqrt2}(e^{9/8}-1)
        >\frac3{2\sqrt2}>1>2/e.
\]

Hence that endpoint test fails. The corrected reference (2) passes by
(8)--(9). The mechanism is exact cancellation of its linear residual,
followed by a full treatment of the remaining cubic and quartic
interactions, including the projected mixed-term cancellation. Merely
stating that a residual has higher algebraic order would not establish
this gain; the positive series estimate (6)--(8) is what closes it.

For the full perturbation class, the initial critical norm need not be
below two, and the initial heat-budget smallness condition need not
hold. No claim is made that the class includes arbitrary initial data,
that further correction steps always succeed, or that this argument
settles the general three-dimensional problem.

The established conclusion is a fully analytic global certificate for
(1), with the complete residual, explicit initial perturbation radius,
and the original NS-B bound (12). Existing certificates are preserved
as independent evidence and are not silently used as premises here.
