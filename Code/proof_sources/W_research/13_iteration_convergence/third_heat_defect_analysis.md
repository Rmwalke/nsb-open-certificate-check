# Exact third heat defect at the fixed datum W

Private local research note. Normalized 2π-periodic three-dimensional torus, viscosity one, real mean-zero divergence-free fields. This calculation concerns the single amplitude a=1. It performs no numerical evolution or time quadrature.

## Construction and signs

Write W=(sin z,sin x,sin y), F=B(W,W), G=B(W,F)+B(F,W), G1=-W/2, G5=G-G1, and H=B(F,F), with B(u,v)=P[(u·∇)v]. Define the nonnegative heat profiles

\[
L_1=\frac{e^{-t}-(1+2t)e^{-3t}}4,\quad
L_5=\frac{e^{-5t}+(2t-1)e^{-3t}}4,\quad
L_6=\frac{(2t^2-2t+1)e^{-4t}-e^{-6t}}4.
\]

Set p=e^{-t}, q=te^{-2t}, d=L1/2, e=L5, f=L6. Then

\[
w_2=(p-d)W-qF+eG_5-fH.
\]

Its full residual r2 has the nine-field expansion

\[
\begin{aligned}
r_2={}&(d^2-2pd)F+qdG+(p-d)eC(W,G_5)-qeC(F,G_5)\\
 &-(p-d)fC(W,H)+qfC(F,H)+e^2B(G_5,G_5)\\
 &+f^2B(H,H)-efC(G_5,H),
\end{aligned}
\]

where C(A,B)=B(A,B)+B(B,A). This is the full projected spatial residual; no modes are discarded.

The next defect and reference are

\[
z_3(t)=\int_0^t e^{(t-s)\Delta}r_2(s)\,ds,
\qquad w_3=w_2-z_3.
\]

The fresh checker constructs all nine spatial fields independently, combines their signed coefficients mode by mode, and convolves each resulting rational exponential polynomial with its exact heat kernel. It finds **342 nonzero Fourier modes**, with maximum squared wave number **22**. Every mode has zero initial value and satisfies its exact forced heat equation; the field is divergence-free. The complete mode profiles are stored in `third_heat_defect_static_result.json`.

## Rigorous global defect bounds

Let X=L²((0,∞); A0), with A0(u)=Σk|ûk|. Exact factorial moments compute every quantity

\[
I_k=\int_0^\infty |\widehat z_{3,k}(t)|^2\,dt\in\mathbb Q_{\ge0}.
\]

Minkowski and integer-square enclosures give

\[
\|z_3\|_X\le\sum_k\sqrt{I_k}
\le\frac{21641807319}{500000000000}
=0.043283614638.
\]

The square integrals retain the cancellations between exponential-polynomial terms inside each mode. The last sum is an upper bound, not an equality for the X norm.

A separate pointwise bound for A0(z3(t)) consists of 31 positive causal field-shell profiles. It retains cancellation inside each scalar residual coefficient. Its complete rational exponential-polynomial representation is under `pointwise_A0_upper_profile` in the JSON; the positivity comes from its causal convolution representation, not from the signs of its expanded exponential coefficients.

## Exact fundamental history and a limitation of the simple global ball

The spatial identities Gfundamental=-W/2 and C(G5,H)fundamental=W/24 imply

\[
(r_2)_{\rm fundamental}=-\left(\frac{qd}2+\frac{ef}{24}\right)W.
\]

Consequently

\[
(z_3)_{\rm fundamental}=-(J_5+J_7)W,
\]

where

\[
J_5=\frac{e^{-t}}{128}[1-(1+2t)e^{-2t}]^2\ge0,
\qquad J_7=\int_0^t e^{-(t-s)}\frac{e(s)f(s)}{24}\,ds\ge0.
\]

The checker verifies the J5 identity algebraically. With J=J5+J7,

\[
\int_0^\infty J^2\,dt
=\frac{64441436014162592029}{34462170976443307130880000}.
\]

Therefore the true X norm also has the lower bound

\[
\|z_3\|_X\ge3\|J\|_{L^2}
\ge\frac{4102348419}{10^{12}}=0.004102348419.
\]

This positive history cannot be reset to zero when moving to a later time interval.

Moreover, L1≤e^{-t}/4 gives

\[
(w_3)_{\rm fundamental}
=[e^{-t}-L_1/2+J_5+J_7]W,
\]

whose scalar coefficient is at least (7/8)e^{-t}. Hence

\[
P_3:=\int_0^\infty A_0(w_3)^2\,dt
\ge\frac{441}{128}>2.
\]

The sufficient simple global contraction-ball hypothesis 2c√P3<1, with c=1/(2√2), thus fails for w3 too. This identifies a limitation of that particular test; it does not imply failure of other continuation criteria or of the underlying solution. The full residual of w3 was not generated or estimated here.

## Block data and verification scope

The JSON key `block_defect_bounds` supplies rigorous bounds on the actual z3 history over the requested grid

0, 1/16, 1/8, 3/16, 1/4, 5/16, 3/8, 7/16, 1/2, 5/8, 3/4, 7/8, 1, 5/4, 3/2, 2, 3, ∞.

These 17 bounds use exact exponential-polynomial antiderivatives, outward rational Taylor/geometric bounds for endpoint exponentials, and modewise L² Minkowski. Their lower bounds use the exact fundamental history. There is no quadrature or stepwise approximation to an evolution equation.

`third_heat_defect_static.py` reports **PASS_THIRD_HEAT_DEFECT_STATIC_CHECKS**, with **2,791 assertions**. This is finite exact arithmetic supporting the written analytic argument. It is not a formal proof-system verification or a numerical simulation of an unknown Navier–Stokes flow. No claim of universal iteration convergence follows.
