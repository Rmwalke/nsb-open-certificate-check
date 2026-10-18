# Enlarged endpoint regularity classes around aW

Private local analytic continuation. No evolution equation or quadrature is computed. Earlier research files remain unchanged.

Let
\[
W=(\sin z,\sin x,\sin y),\qquad
\chi(h)=A_{-1}(h)=\sum_{k\ne0}\frac{|\widehat h_k|}{|k|},
\]
on the normalized \(2\pi\)-periodic three-dimensional torus with viscosity one. The initial datum is \(u_0=aW+h\), where \(a\) is real and \(h\) is smooth, real, mean-zero, and divergence-free, with arbitrary Fourier support.

The following closed classes are now certified to have a global smooth solution and a cutoff-independent NS-B bound:
\[
\boxed{|a|\le1,\quad \chi(h)\le1/10,}                         \tag{1A}
\]
\[
\boxed{|a|\le11/10,\quad \chi(h)\le1/40,}                     \tag{1B}
\]
\[
\boxed{|a|\le8/7,\quad \chi(h)\le1/100.}                      \tag{1C}
\]
Two contained classes have smaller error barriers and useful stronger bounds:
\[
|a|\le1,\quad\chi(h)\le2/25,\qquad
|a|\le9/8,\quad\chi(h)\le1/100.                              \tag{1D}
\]

These genuinely enlarge the previously established union
\[
\{|a|\le1,\ \chi(h)\le1/20\}
\ \cup\
\{|a|\le11/10,\ \chi(h)\le1/100\}.
\]
They are endpoint existence and NS-B certificates. No ordinary Picard convergence claim is inferred for these enlarged classes.

## 1. What changes in the endpoint test

The same exact second corrected reference \(w_{2,a}\) is used; no new evolution or larger Fourier reference is needed. Three pieces of information are retained more accurately:

1. The initial transformed error is \(\Psi(\chi)=\chi e^{-\chi/2}\), not its coarser upper bound \(\chi\).
2. Its multiplier uses the exact second-reference budget \(P_2=\int A_0(w_{2,a})^2\), enclosed by rational arithmetic, rather than \(P_1\).
3. The endpoint transform is allowed its full increasing range below two, instead of imposing an error cap of one.

Write
\[
b=A_{-1}(u-w_{2,a}),\quad h_e=A_1(u-w_{2,a}),\quad
U=A_0(w_{2,a}),\quad R\ge A_{-1}(r_2).
\]
The established complete error inequality is
\[
b'+(1-b/2)h_e\le U\sqrt{b h_e}+R.                            \tag{2}
\]
All residual modes and the complete initial perturbation occur here.

On \(0\le b<2\), optimize the right side minus its \(h_e\) term and set
\(\Psi(b)=be^{-b/2}\). Since
\(\Psi'(b)=e^{-b/2}(1-b/2)\in(0,1]\),
\[
(\Psi(b))'\le \tfrac14 U^2\Psi(b)+R.
\]
For \(P=\int_0^\infty U^2\), a sufficient gate with any \(0\le r<B<2\) is
\[
\boxed{
M:=e^{P/4}\Psi(r)+
\int_0^\infty
 e^{\frac14\int_t^\infty U(s)^2ds}R(t)\,dt
<\Psi(B),\qquad b(0)\le r.
}                                                               \tag{3}
\]
The inverse of \(\Psi\) is used only on \([0,2)\). Its first-contact argument gives \(b(t)<B\) globally. The argument first applies to finite Galerkin systems containing the reference support; it does not assume global regularity beforehand.

The two terms of (3) may be bounded differently. Below, the initial term uses exact \(P_2\), whereas the residual weight uses the older valid remaining-time bound through \(U_1\). There is no assumption that these different upper bounds are equal.

## 2. Exact reference amplitude and interval monotonicity

The second reference is
\[
w_{2,a}=(ae^{-t}-a^3L_1/2)W-a^2te^{-2t}F
                         +a^3L_5C-a^4L_6H,                 \tag{4}
\]
where
\[
F=B(W,W),\quad
C=\tfrac12(\sin z\cos2y,\sin x\cos2z,\sin y\cos2x),\quad
H=B(F,F),
\]
and
\[
\begin{split}
L_1&=[e^{-t}-(1+2t)e^{-3t}]/4,\\
L_5&=[e^{-5t}+(2t-1)e^{-3t}]/4,\\
L_6&=[(2t^2-2t+1)e^{-4t}-e^{-6t}]/4.
\end{split}
\]
These are nonnegative heat-convolution profiles with zero initial value. The reference has support \(|k|^2\le6\) and starts exactly from \(aW\), so its entire initial error is \(h\).

Let \(m=|a|\). The disjoint reference shells and their known Fourier norms give
\[
U_2=3me^{-t}+3m^2te^{-2t}
 -\tfrac32m^3(L_1-L_5)+\sqrt3m^4L_6.                         \tag{5}
\]
This identity holds for both signs of \(a\) when \(m^2\le8\), because the fundamental coefficient keeps the sign of \(a\).

The established profile inequalities are
\[
0\le L_1-L_5\le L_1\le e^{-t}/4,\qquad L_1-L_5\ge2L_6.
\]
The second follows by comparing the ordered decay rates in
\[
L_1-L_5=4k_1*k_3*k_3*k_5,\qquad
2L_6=4k_4*k_4*k_4*k_6,
\quad k_\lambda(t)=e^{-\lambda t}\mathbf1_{t\ge0}.
\]
Consequently,
\[
U_2\le U_1:=3me^{-t}+3m^2te^{-2t}\qquad(m^2\le3).             \tag{6}
\]

For the use of an exact endpoint budget over a whole amplitude interval, monotonicity is needed separately. Differentiating (5) gives
\[
\partial_m U_2
\ge (3-\tfrac98m^2)e^{-t}\ge0
\qquad(m^2\le8/3).                                           \tag{7}
\]
The omitted derivative terms are nonnegative. The same proof holds after replacing \(\sqrt3\) by a rational upper enclosure, since the extra term multiplies \(L_6\ge0\).

Every declared candidate endpoint is at most \(5/4\), which lies strictly within both monotonicity ranges. Thus a rational upper enclosure \(P_2^+(A)\) computed at amplitude \(A\) bounds the exact \(P_2(m)\) for all \(0\le m\le A\).

The new checker encloses \(\sqrt3\) between adjacent multiples of \(10^{-18}\), forms the square of the upper version of (5), and evaluates its exact exponential-polynomial moment. The numerical orientations relevant here are
\[
\begin{array}{c|c}
 A & P_2^+(A)\\ \hline
1 & 6.569765510001933\ldots\\
11/10 & 8.207493783861516\ldots\\
9/8 & 8.652517455190253\ldots\\
8/7 & 8.979357269363748\ldots
\end{array}
\]
The rational values, rather than these decimal displays, are used in all gates.

## 3. Complete residual with the useful signed groups retained

For \(m\ge0\), let
\[
x=me^{-t},\quad y=m^2te^{-2t},\quad
p=m^3L_1/2,\quad q=m^3L_5,\quad r_6=m^4L_6.
\]
With \(K_{AB}=B(A,B)+B(B,A)\) and
\(G=K_{WF}=-W/2+C\), the full residual is
\[
\begin{split}
r_2={}&(p^2-2xp)F+ypG+(x-p)qK_{WC}-(x-p)r_6K_{WH}\\
 &-yqK_{FC}+yr_6K_{FH}+q^2B(C,C)-qr_6K_{CH}+r_6^2B(H,H).
\end{split}                                                   \tag{8}
\]
For negative amplitudes, taking the grouped coefficient magnitudes produces the same scalar functions of \(m=|a|\).

The previously verified complete Fourier caps are reused as mathematical input:

| Field | \(A_{-1}\) upper cap |
|---|---:|
| \(F\) | \(17/8\) |
| \(K_{WC}\) | \(171/100\) |
| \(K_{WH}\) | \(27/20\) |
| \(G\) | \(87/40\) |
| \(K_{FC}\) | \(6/5\) |
| \(K_{FH}\) | \(183/100\) |
| \(B(C,C)\) | \(11/25\) |
| \(K_{CH}\) | \(3/4\) |
| \(B(H,H)\) | \(9/25\) |

With these caps denoted \(c_j\), define
\[
R_m=\sum_{j=1}^9 c_j f_j,\quad
(f_1,\ldots,f_9)=
(2xp-p^2,(x-p)q,(x-p)r_6,yp,yq,yr_6,q^2,qr_6,r_6^2).          \tag{9}
\]
Then \(R_m\ge A_{-1}(r_2)\). All nine scalar profiles are nonnegative in the candidate range. The subtractive terms inside the first three profiles are preserved before taking their Fourier-field norm caps.

The first three profiles have amplitude form \(m^4A-m^6D\), \(m^4A-m^6D\), and \(m^5A-m^7D\), with \(D/A\le1/16,1/8,1/8\), respectively. Their derivative lower factors are
\(4-6m^2/16\), \(4-6m^2/8\), and \(5-7m^2/8\), all positive for \(m^2\le3\). The other profiles are nonnegative pure powers. Hence \(R_m\) is pointwise nondecreasing throughout the required amplitude intervals.

## 4. Exact moments and the positive residual-weight series

Although an expanded \(R_m\) has some negative exponential-polynomial coefficients, its function is nonnegative by (9). For each nonnegative integer \(n\),
\[
I_n(m)=\int_0^\infty e^{-2nt}R_m(t)dt
\]
is therefore nonnegative and nonincreasing in \(n\). Each value is an exact rational combination of
\[
\int_0^\infty t^d e^{-(\lambda+2n)t}dt
=\frac{d!}{(\lambda+2n)^{d+1}}.                              \tag{10}
\]

The remaining reference budget has the valid bound
\[
\exp\!\left(\frac14\int_t^\infty U_2^2\right)
\le e^{\Delta_m}\exp(c_m e^{-2t}),
\quad
c_m=\frac98m^2,\quad
\Delta_m=\frac12m^3+\frac9{128}m^4.                           \tag{11}
\]
This follows from \(U_2\le U_1\) and the lower heat component of the past \(U_1^2\) integral. It does not claim the false general lower bound \(U_2\ge3me^{-t}\).

Expanding the final positive exponential and using Tonelli gives
\[
\int_0^\infty e^{\frac14\int_t^\infty U_2^2}R_m(t)dt
\le e^{\Delta_m}\sum_{n=0}^\infty\frac{c_m^n}{n!}I_n(m).
\]
Set \(t_n=c_m^nI_n/n!\). Since \(I_{n+1}\le I_n\), the fresh bound is
\[
\boxed{
M_{\rm res}^+(m)=e^{\Delta_m}
\left[\sum_{n=0}^{24}t_n+\frac{t_{25}}{1-c_m/26}\right].
}                                                               \tag{12}
\]
All displayed exponentials are replaced by positive rational Taylor-series upper bounds with explicit geometric tails.

For a class endpoint \((A,r)\), the resulting gate is
\[
\boxed{
M^+(A,r)=
r\exp\!\left(P_2^+(A)/4-r/2\right)+M_{\rm res}^+(A)
<\Psi(B).
}                                                               \tag{13}
\]
This retains the exact initial transform \(\Psi(r)\). The right side is certified from below by
\(B/\overline{e^{B/2}}\). The inequalities are strict at the closed class endpoints, so the claimed non-strict input bounds are valid.

## 5. Passing endpoint classes

The following numbers are decimal orientations. Each row has a strictly positive rational margin in the result JSON.

| \(A\) | Initial radius \(r\) | Error cap \(B\) | \(M^+(A,r)\) | Certified margin below \(\Psi(B)\) |
|---:|---:|---:|---:|---:|
| 1 | \(2/25\) | \(6/5\) | 0.636555727041 | \(>0.02201823627\) |
| 1 | \(1/10\) | \(9/5\) | 0.730918119157 | \(>0.00090726837\) |
| \(11/10\) | \(1/40\) | \(6/5\) | 0.656895071385 | \(>0.00167889192\) |
| \(9/8\) | \(1/100\) | \(6/5\) | 0.636557679768 | \(>0.02201628354\) |
| \(8/7\) | \(1/100\) | \(8/5\) | 0.714761998246 | \(>0.00416434434\) |

The initial radius \(1/10\) around the unit segment is the tightest displayed comparison. Its margin is still an exact positive rational, rather than a rounded numerical observation.

The signed residual totals needed for continuation are bounded by
\[
\begin{array}{c|c}
 A & \int_0^\infty R_A(t)\,dt\\ \hline
1 & 0.100072634381284\ldots\\
11/10 & 0.148566222682886\ldots\\
9/8 & 0.163091475625743\ldots\\
8/7 & 0.174113083988406\ldots
\end{array}
\]
Again, the result JSON stores their exact rational values.

## 6. Explicit original NS-B consequence

The sharper continuation estimate proved independently in *enlarged_budget_consequence.md* retains the remaining reference energy in the endpoint envelope. Here is its analytic form, so the NS-B consequence is explicit even without rounding to an integer exponent.

Define, for \(0\le s<2\),
\[
F_*(s)=-2s-4\log(1-s/2),\qquad
\delta_B=1-B/2,\qquad
C=F_*(r)+\frac{B}{\delta_B}R_*,
\]
where \(R_*=\int R_A\). For any passing row let
\[
\boxed{
G=\left[\sqrt{P_2^+(A)}+\sqrt{F_*(B)}
                  +\sqrt{F_*(B)+C}\right]^2.
}                                                               \tag{14}
\]
Then every Galerkin cutoff containing the second reference satisfies
\[
\int_0^\infty A_0(u_K)^2dt\le G.
\]

For completeness, (3) contains the stronger time-dependent information
\[
\Psi(b(t))\le
\Psi(B)\exp\!\left(-\frac14\int_t^\infty U_2^2\right).
\]
Define \(\beta(t)\in[0,B]\) by equality. Then
\(\beta'=U_2^2\beta/[4(1-\beta/2)]\), and
\[
\int_0^\infty U_2^2\frac{b^2}{(1-b/2)^2}dt\le4F_*(B).
\]
Multiplying (2) by \(b/(1-b/2)=F_*'(b)\), integrating first over a finite interval, and applying Cauchy–Schwarz yields
\[
Q_T:=\int_0^T b h_e
 \le C+\sqrt{4F_*(B)Q_T}.
\]
Solving before letting \(T\) increase gives
\(\sqrt Q\le\sqrt{F_*(B)}+\sqrt{F_*(B)+C}\). Since
\(A_0(u-w_2)^2\le b h_e\), Minkowski proves (14). No finiteness of \(Q\) is assumed in advance.

The reference support is \(|k|\le\sqrt6\). For \(K\ge\sqrt6\), projection contracts the complete residual and initial critical error, so the entire comparison is uniform in \(K\). The coefficient-one \(H^1\) inequality and the energy identity therefore give
\[
\int_0^\infty\|\nabla u_K(t)\|_2^4dt
 \le\frac{E_0X_0}{2}e^G,\qquad
E_0=\|u_0\|_2^2,\quad X_0=\|\nabla u_0\|_2^2.                 \tag{15}
\]

For \(K<\sqrt6\), use energy directly:
\[
\sup_t\|\nabla u_K\|_2^2\le6E_0,\qquad
\int_0^\infty\|\nabla u_K\|_2^2dt\le E_0/2.
\]
Hence
\[
\int_0^\infty\|\nabla u_K\|_2^4dt
 \le3E_0^2\le3E_0X_0.
\]
Every reported budget can be rounded upward to an exponent \(L\) with \(e^L\ge6\), giving
\[
\boxed{\sup_K\int_0^\infty\|\nabla u_K(t)\|_2^4dt
       \le\frac{E_0X_0}{2}e^L,\qquad L\ge G.}                 \tag{16}
\]
The lower cutoffs are not falsely claimed to satisfy the same reference-amplitude budget. The full initial energy and higher Sobolev norms need not be small in these critical neighborhoods. The uniform estimates give a global strong limit and smooth continuation for smooth initial data.

The fresh rational logarithm-and-square-root check in *enlarged_budget_static.py* certifies the following concrete choices in (16):

| Amplitude bound \(A\) | Critical radius \(r\) | Integer NS-B exponent \(L\) |
|---:|---:|---:|
| 1 | \(2/25\) | 25 |
| 1 | \(1/10\) | 59 |
| \(11/10\) | \(1/40\) | 29 |
| \(9/8\) | \(1/100\) | 30 |
| \(8/7\) | \(1/100\) | 49 |

Its exact enclosures are stored in *enlarged_budget_static_result.json*, with the endpoint source hash checked before use. The retained-history estimate also reduces the common exponent for the two previously proved classes from 81 to 24; that consequence retains their original error cap one and does not assume an unproved combined amplitude-radius rectangle.

## 7. Strict enlargement rather than a different decomposition

The new union contains both old classes: (1A) contains the old unit segment tube, and (1C) contains the old \(11/10\) segment tube.

A concrete datum outside both old classes is
\[
u_0=(8/7)W.
\]
For any old center \(aW\) with \(|a|\le1\), its required perturbation has
\[
\chi((8/7-a)W)\ge3/7>1/20,
\]
because \(A_{-1}(W)=3\). For any old center with \(|a|\le11/10\),
\[
\chi((8/7-a)W)\ge3(8/7-11/10)=9/70>1/100.
\]
Thus no alternate allowed old amplitude decomposition places this datum in either old tube. Class (1C) includes it with zero perturbation.

The increased unit perturbation radius also adds transverse data, not merely larger multiples of \(W\). For example
\[
u_0=W+\frac15\sin(2x)e_2
\]
has \(\chi(h)=1/10\) and disjoint Fourier support between \(h\) and \(W\). For any scalar old center \(aW\),
\[
\chi(u_0-aW)=3|1-a|+1/10.
\]
This is greater than either old radius for every allowed old amplitude, while (1A) contains the datum.

These comparisons concern the stated earlier analytic classes. They do not assert that every new point was unknown to every possible independent certificate.

## 8. Bounded verification and limits

The fresh standard-library program *enlarged_endpoint_gate.py* constructs the exponential polynomials anew and uses exact rational factorial moments. It imports no historical research module and executes no historical checker. The previously established complete Fourier-field caps are reused as proved mathematical input, with source hashes recorded.

The finite declared candidates were
\[
(1,2/25),(1,1/10),(1,11/100),
(11/10,1/40),(11/10,1/20),
(9/8,1/100),(9/8,1/40),
(8/7,1/100),(7/6,1/100),(6/5,1/100),(5/4,0).
\]
The single \(8/7\) candidate was a bounded follow-up suggested before finalization; no threshold search followed. The five passing rows are listed above. The other cases were not certified by their computed positive upper majorants. Such noncertification proves neither a singularity nor failure of the actual endpoint gate or a different reference.

The checker reports **PASS_ENLARGED_ENDPOINT_REGULARITY_CLASSES**, with **465 exact assertions**. Its frozen output is *enlarged_endpoint_gate_result.json*, SHA-256
\[
\texttt{661a9c1c426b59e5f9d0a9477628f03f3b96a1ec8f86dc69f4c7d5e919ae4f0f}.
\]

The work uses no numerical evolution, quadrature, web access, upload, or modification of earlier-stage artifacts. It is an analytic certificate supported by finite exact arithmetic, not a formal proof-assistant verification. The arbitrary-data problem and convergence of ordinary Picard iteration on these larger tubes remain outside this result.
