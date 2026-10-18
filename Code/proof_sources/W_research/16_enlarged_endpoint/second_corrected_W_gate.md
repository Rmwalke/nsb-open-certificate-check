# Second residual correction: larger amplitude and critical neighborhoods

Private analytic result, 2026-09-30 (Asia/Kolkata). This note uses new exact
Fourier coefficient checks and a new standalone calculation of rational
factorial moments and positive-series tails. No historical program was
executed, no ODE/PDE was numerically integrated, no external source was
used, and earlier files remain unchanged.

The second residual correction proves either of the following sufficient
conditions for smooth real, mean-zero, divergence-free initial data:

\[
 \boxed{u_0=aW+h,\quad |a|\le11/10,\quad A_{-1}(h)\le1/100,}
                                                               \tag{1A}
\]
\[
 \boxed{u_0=aW+h,\quad |a|\le1,\quad A_{-1}(h)\le1/20,}
 \qquad W=(\sin z,\sin x,\sin y).                \tag{1B}
\]

The perturbation may have arbitrary full Fourier support. Both classes
have a global smooth solution and the explicit original NS-B bound

\[
 \boxed{\sup_K\int_0^\infty\|\nabla u_K(t)\|_2^4dt
                 \le\frac{E_0X_0}{2}e^{81},
 \quad E_0=\|u_0\|_2^2,\quad X_0=\|\nabla u_0\|_2^2.}        \tag{2}
\]

This improves the previous analytically certified amplitude segment and
critical perturbation radii. It does not assert arbitrary-data regularity
or optimal constants.

## 1. Exact second reference and its time profiles

Use viscosity one, normalized measure on the \(2\pi\)-torus,
\(A=-\Delta\), and \(B(u,v)=\mathbb P((u\cdot\nabla)v)\).
The previous exact Fourier fields are

\[
 F=B(W,W),\quad G=B(W,F)+B(F,W)=-W/2+C,\quad
 C=\tfrac12(\sin z\cos2y,\sin x\cos2z,\sin y\cos2x),
 \quad H=B(F,F).
\]

Their squared frequency shells are respectively \(2\), \(1,5\),
\(5\), and \(6\). The first corrected reference is
\(w_1=a e^{-t}W-a^2t e^{-2t}F\), with full residual
\(r_1=-a^3t e^{-3t}G+a^4t^2e^{-4t}H\).

Define positive convolution profiles

\[
 L_1=\frac{e^{-t}-(1+2t)e^{-3t}}4,\quad
 L_5=\frac{e^{-5t}+(2t-1)e^{-3t}}4,
\]
\[
 L_6=\frac{(2t^2-2t+1)e^{-4t}-e^{-6t}}4.           \tag{3}
\]

They vanish initially and solve
\((\partial_t+1)L_1=(\partial_t+5)L_5=t e^{-3t}\) and
\((\partial_t+6)L_6=t^2e^{-4t}\). Positivity follows from their
causal heat convolutions, despite the minus signs in the displayed
exponential-polynomial expressions.

The second full linear residual response is

\[
 z_2=(a^3/2)L_1W-a^3L_5C+a^4L_6H,
 \qquad z_2'+Az_2=r_1,\quad z_2(0)=0.
\]

Thus the corrected reference is exactly

\[
 \boxed{w_2=(a e^{-t}-a^3L_1/2)W-a^2t e^{-2t}F
                     +a^3L_5C-a^4L_6H.}           \tag{4}
\]

It has support inside \(|k|\le\sqrt6\) and initial value \(aW\).
Its complete residual is
\(r_2=-B(w_1,z_2)-B(z_2,w_1)+B(z_2,z_2)\).

## 2. The reference amplitude does not grow under this correction

Write \(m=|a|\). Since \(0\le L_1\le e^{-t}/4\), the fundamental
coefficient in (4) retains its sign for \(m\le\sqrt8\). The four
shell supports are disjoint, and
\(A_0(W)=A_0(F)=3\), \(A_0(C)=3/2\), \(A_0(H)=\sqrt3\).
Consequently the exact reference norm is

\[
 U_2=A_0(w_2)=3m e^{-t}+3m^2t e^{-2t}
               -\tfrac32m^3(L_1-L_5)+\sqrt3m^4L_6.            \tag{5}
\]

The negative fundamental correction has been preserved in (5); replacing
all four coefficients by a triangle bound would discard a useful
cancellation.

To compare (5) with \(U_1=3m e^{-t}+3m^2t e^{-2t}\), let
\(k_\lambda(t)=e^{-\lambda t}\mathbf1_{t\ge0}\). Direct resolvent
identities give

\[
 L_1-L_5=4k_1*k_3*k_3*k_5,\qquad
 2L_6=4k_4*k_4*k_4*k_6.
\]

Each decay rate in the first ordered list \((1,3,3,5)\) is no larger
than the corresponding rate in \((4,4,4,6)\). Positivity of convolution
therefore proves \(L_1-L_5\ge2L_6\). It follows that

\[
 U_2-U_1\le m^3(\sqrt3m-3)L_6\le0
                  \quad\text{for }0\le m\le\sqrt3.           \tag{6}
\]

In particular, with \(P_i=\int_0^\infty U_i^2dt\),

\[
 P_2\le P_1(m)=\frac92m^2+2m^3+\frac9{32}m^4.      \tag{7}
\]

More importantly, (6) also compares every remaining time budget:
\(\int_t^\infty U_2^2\le\int_t^\infty U_1^2\). Thus the endpoint
weight satisfies

\[
 \boxed{\exp\left(\tfrac14\int_t^\infty U_2^2\right)
 \le e^{\Delta_m}\exp(c_m e^{-2t}),\quad
 c_m=9m^2/8,\quad \Delta_m=m^3/2+9m^4/128.}         \tag{8}
\]

This uses \(U_1\ge3m e^{-t}\) to subtract the heat part of its past
budget. It does not assert the generally unjustified lower inequality
\(U_2\ge3m e^{-t}\).

## 3. A complete positive residual majorant

For nonnegative amplitude, introduce nonnegative scalar functions

\[
 x=m e^{-t},\quad y=m^2t e^{-2t},\quad
 p=m^3L_1/2,\quad q=m^3L_5,\quad r=m^4L_6.
\]

Then \(w_1=xW-yF\) and \(z_2=pW-qC+rH\). Write
\(K_{AB}=B(A,B)+B(B,A)\). Exact expansion gives

\[
\begin{aligned}
 r_2={}&(p^2-2xp)F+ypG+(x-p)qK_{WC}-(x-p)rK_{WH}\\
       &-yqK_{FC}+yrK_{FH}+q^2B(C,C)-qrK_{CH}+r^2B(H,H).
\end{aligned}                                                   \tag{9}
\]

For signed amplitudes, taking the coefficient magnitudes produces the
same functions of \(m=|a|\). Since \(p/x\le m^2/8\), every scalar
coefficient magnitude used below is nonnegative for \(m\le\sqrt3\).

Fresh exact Fourier enumeration, including Leray projection and every
output mode, proves the following convenient rational norm caps:

| Field | Strict upper bound for \(A_{-1}\) | Nonzero modes |
|---|---:|---:|
| \(F\) | \(17/8\) | 12 |
| \(K_{WC}\) | \(171/100\) | 48 |
| \(K_{WH}\) | \(27/20\) | 72 |
| \(G\) | \(87/40\) | 18 |
| \(K_{FC}\) | \(6/5\) | 60 |
| \(K_{FH}\) | \(183/100\) | 120 |
| \(B(C,C)\) | \(11/25\) | 48 |
| \(K_{CH}\) | \(3/4\) | 162 |
| \(B(H,H)\) | \(9/25\) | 120 |

Denote these nine caps in table order by \(c_j\), and define

\[
 (f_1,\ldots,f_9)=
 (2xp-p^2,(x-p)q,(x-p)r,yp,yq,yr,q^2,qr,r^2),
\]
\[
 \boxed{R_m(t):=\sum_{j=1}^9 c_j f_j(t)
                  \ge A_{-1}(r_2(t)).}            \tag{10}
\]

The negative terms inside the first three scalar profiles are retained
before taking the norm bounds. Other cross-field cancellations are not
needed for the proof.

The majorant is pointwise nondecreasing in \(m\) for
\(0\le m\le\sqrt3\). To check the only nontrivial cases, their forms
are \(m^4A-m^6D\), \(m^4A-m^6D\), and \(m^5A-m^7D\), with
\(D/A\le1/16,1/8,1/8\), respectively. Their derivatives have positive
factors bounded below by
\(4-6m^2/16\), \(4-6m^2/8\), and \(5-7m^2/8\).
The other profiles are nonnegative pure powers of \(m\). This proves
the entire amplitude intervals below, rather than only their endpoints.

## 4. Exact factorial moments and positive-series gate

For a fixed rational \(m\), every \(f_j\), and hence \(R_m\), is a
finite exponential polynomial with rational coefficients. Its weighted
moments are computed exactly by

\[
 I_n(m)=\int_0^\infty e^{-2nt}R_m(t)dt,
 \qquad
 \int_0^\infty t^d e^{-(\lambda+2n)t}dt
                  =\frac{d!}{(\lambda+2n)^{d+1}}.              \tag{11}
\]

Although individual exponential-polynomial coefficients can be negative,
the function \(R_m\) is nonnegative by its construction in (10).
Therefore \(I_n\ge0\) and \(I_{n+1}\le I_n\).

Combining (8), (10), and the positive exponential series gives

\[
 M_{\rm res}(m)\le e^{\Delta_m}S_m,\qquad
 S_m=\sum_{n=0}^\infty\frac{c_m^n}{n!}I_n(m).       \tag{12}
\]

Tonelli applies to the nonnegative integrands. The factorial identities
in (11) are exact integrals of specified functions, not numerical ODE
or PDE integration.

Let \(t_n=c_m^nI_n/n!\). Then
\(t_{n+1}/t_n\le c_m/(n+1)\), and the finite bound used here is

\[
 S_m\le\sum_{n=0}^{8}t_n+\frac{t_9}{1-c_m/10}.     \tag{13}
\]

The new standalone scalar checker verifies the following exact rational
comparisons using (11)--(13):

| Quantity | \(m=11/10\) | \(m=1\) |
|---|---:|---:|
| \(\int R_m\) | \(<3/20\) | \(<101/1000\) |
| \(S_m\) | \(<27/125\) | \(<271/2000\) |
| \(e^{\Delta_m}\) | \(<54/25\) | \(<177/100\) |
| \(M_{\rm res}(m)\) | \(<47/100\) | \(<6/25\) |
| \(e^{P_1(m)/4}\) | \(<17/2\) | \(<11/2\) |

All exponential bounds use positive Taylor series with geometric tails.
For example, for the arguments here the checker uses

\[
 e^x\le\sum_{n=0}^{16}\frac{x^n}{n!}
        +\frac{x^{17}/17!}{1-x/18},\qquad 0\le x<18.
\]

## 5. Both endpoint gates pass with the same error barrier

For \(u_0=aW+h\), the initial error of (4) is exactly \(h\). Put
\(b_0=A_{-1}(h)\) and \(\Psi(b)=be^{-b/2}\). The actual endpoint
majorant is bounded by

\[
 M\le e^{P_1(m)/4}b_0+M_{\rm res}(m).
\]

In class (1A), the preceding table gives

\[
 M<\frac{47}{100}+\frac{17/2}{100}
             =\frac{111}{200}<\frac35.
\]

In class (1B), it gives

\[
 M<\frac6{25}+\frac{11/2}{20}
             =\frac{103}{200}<\frac35.
\]

The same Taylor-tail method verifies \(e^{1/2}<5/3\), so
\(\Psi(1)=e^{-1/2}>3/5\). The first-contact endpoint proof therefore
gives the uniform complete-error barrier

\[
 \boxed{A_{-1}(u-w_2)<1.}                          \tag{14}
\]

This applies first to every Galerkin cutoff containing the reference
support. Projecting the full residual can only reduce its critical
norm, and projecting the initial perturbation can only reduce \(b_0\).
Thus no global full-PDE solution is assumed in deriving (14).

## 6. Explicit continuation and the NS-B constant 81

Both classes have

\[
 P_2\le P_1(11/10)=2726009/320000<213/25<9.
\]

The sum of initial error and total residual is below \(4/25\):
in (1A) use \(1/100+3/20=4/25\); in (1B) use
\(1/20+101/1000<4/25\). Let \(h_e=A_1(u-w_2)\).
Returning to the unsplit error inequality with \(b<1\) gives

\[
 b'+\tfrac12h_e\le U_2\sqrt{h_e}+R_m
             \le\tfrac14h_e+U_2^2+R_m.
\]

Consequently

\[
 \int_0^\infty h_e dt
 \le4P_2+4\left(b_0+\int_0^\infty R_mdt\right)
 <\frac{868}{25}<35.                              \tag{15}
\]

Since \(A_0(u-w_2)^2\le b h_e\), the error squared-amplitude integral
is below 35. Minkowski in time proves, for every higher cutoff,

\[
 \int_0^\infty A_0(u_K)^2dt
       <(\sqrt9+\sqrt{36})^2=81.                  \tag{16}
\]

Lower cutoffs \(K<\sqrt6\) can be covered separately without claiming
that the full perturbation energy is small. In (1A), Fourier triangle
inequalities give

\[
 \|P_Ku_0\|_2
 \le(11/10)\sqrt{3/2}+\sqrt6/100
       <(11/10)(5/4)+(5/2)/100=7/5.
\]

In (1B), the corresponding bound
\(\sqrt{3/2}+\sqrt6/20<5/4+(5/2)/20<7/5\) also holds.
The exact lattice enumeration through squared frequency six gives

\[
 S_2(\sqrt6)=6+12/2+8/3+6/4+24/5+24/6=749/30<25.
\]

Thus the low-cutoff energy estimate yields
\(\int A_0(u_K)^2\le S_2(K)\|P_Ku_0\|_2^2/2<49/2<81\).
Combining this with (16) covers every cutoff. Testing the equation at
\(H^1\), and using the ordinary energy equality, gives

\[
 \sup_{K,t}\|\nabla u_K(t)\|_2^2\le X_0e^{81},\qquad
 \sup_K\int_0^\infty\|\nabla u_K(t)\|_2^4dt
                     \le\frac{E_0X_0}{2}e^{81}.
\]

The uniform bounds give a global strong limit, and smooth initial data
give global smoothness. The full datum's \(E_0,X_0\) can be arbitrarily
large within these critical perturbation neighborhoods; only the
low-frequency projections were bounded in the preceding paragraph.

## 7. Verification, comparison, and limitation

The new files `second_correction_fourier_check.py` and its result JSON
verify the finite Fourier fields, projections, profile equations, and
rational norm caps. The new `second_corrected_gate_static.py` and its
result JSON verify 69 exact scalar assertions, including the factorial
moments, positive-series remainder, both radii, and exponent 81.
Neither checker executes or imports an earlier research program.

Relative to the previous first-correction certificate, the amplitude
endpoint increases from one to \(11/10\), its accompanying critical
radius increases from \(1/500\) to \(1/100\), and the radius around
the unit-amplitude segment increases to \(1/20\). The explicit analytic
NS-B exponent decreases from 6541 to 81. The older numerical W
certificate retains its independent, substantially sharper bound for
its own certified neighborhood.

These are concrete gains from the second correction, with complete
residual accounting. There is no claim that repeating this correction
must converge for arbitrary data, or that failure of one of these
sufficient gates indicates a singularity. The general three-dimensional
question remains open.
