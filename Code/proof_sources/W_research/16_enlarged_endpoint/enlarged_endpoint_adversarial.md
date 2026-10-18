# Independent audit: enlarged endpoint regularity classes

This note audits new classes around \(aW\), where
\(W=(\sin z,\sin x,\sin y)\), using the previously constructed second
corrected reference \(w_{2,a}\). The old classes were
\(|a|\le1,\ \chi(h)\le1/20\) and
\(|a|\le11/10,\ \chi(h)\le1/100\), with error cap one.
Here \(\chi=A_{-1}\), the initial datum is \(aW+h\), and perturbations are
smooth, real, mean-zero and divergence-free. The reference still begins
at \(aW\); the complete \(h\) remains the initial error.

The audit passes. In particular, the critical radius at \(|a|\le1\)
increases to \(1/10\), and the amplitude endpoint at radius \(1/100\)
increases to \(8/7\). These are strict enlargements of the union of the
two named earlier classes, as verified below. No previous source file was
executed, imported or modified for this audit.

## 1. Monotone branch and the initial critical radius

The error inequality is

\[
b'+\delta(b)h_e\le U\sqrt{b h_e}+R,
\qquad \delta(b)=1-b/2,
\]

where \(b=A_{-1}(u-w_2)\), \(h_e=A_1(u-w_2)\),
\(U=A_0(w_2)\), and \(R\) bounds the full critical residual.
Write \(A(t)=\int_0^tU^2\), \(P=A(\infty)\), and
\(\Psi(b)=be^{-b/2}\). The endpoint transform is increasing precisely
on \([0,2]\); all inversion and comparison here use the branch \([0,2)\).

For \(0\le r<B<2\), a sufficient certificate is

\[
e^{P/4}\Psi(r)+
\int_0^\infty e^{(P-A(t))/4}R(t)\,dt<\Psi(B).
\tag{1}
\]

The actual \(b(0)\) is at most \(r\), so
\(\Psi(b(0))\le\Psi(r)\). A first-contact argument at \(B\) gives
\(b(t)<B\), preventing approach to the turning point two. Replacing
\(\Psi(r)\) by \(r\) would also be safe but unnecessarily weaker.
Applying this monotonicity to an initial radius at or beyond two would
not be justified. No inverse on the decreasing branch is used.

The two terms in (1) may be bounded by different valid upper estimates.
In particular, exact \(P_2\) can be used for the initial term while the
residual term uses the previously proved \(U_2\le U_1\) tail bound.
This does not require a lower bound on \(A_2(t)\), nor equality between
the two upper estimates.

## 2. Signed amplitudes and safe monotonicity ranges

Let \(m=|a|\). The exact reference amplitude is

\[
U_2=3me^{-t}+3m^2te^{-2t}
-\frac32m^3(L_1-L_5)+\sqrt3m^4L_6,
\]

where the nonnegative heat profiles satisfy
\(L_1\le e^{-t}/4\), \(L_1-L_5\ge2L_6\ge0\).
The fundamental reference coefficient retains the sign of \(a\) when
\(m^2\le8\). Thus the displayed expression is valid for both signs of
\(a\) throughout every range considered here.

The previous domination \(U_2\le U_1=3me^{-t}+3m^2te^{-2t}\) is valid
for \(m\le\sqrt3\). Moreover

\[
\partial_m U_2\ge
\left(3-\frac98m^2\right)e^{-t}\ge0
\quad\text{when }m^2\le8/3.
\]

Indeed \(0\le L_1-L_5\le L_1\le e^{-t}/4\), and the remaining derivative
terms are nonnegative. This proof also applies when \(\sqrt3\) is replaced
by a rational upper bound: its contribution multiplies \(L_6\ge0\).
Consequently the exact \(P_2\) upper enclosure is monotone in this range.
Candidates through \(m=5/4\) lie strictly within it.

The complete nine-field residual majorant retains the three grouped
positive differences rather than bounding their summands separately.
Their amplitude structures are
\(m^4A-m^6B\), \(m^4A-m^6B\), and \(m^5A-m^7B\), with coefficient
ratios \(B/A\le1/16,1/8,1/8\), respectively. Differentiation gives lower
factors \(4-6m^2/16\), \(4-6m^2/8\), and \(5-7m^2/8\), all positive for
\(m^2\le3\). The other profiles are nonnegative pure powers. Their
magnitudes are the same for signed \(a\). Thus the full residual majorant
is pointwise nondecreasing over the required amplitude intervals, even
though an expanded exponential-polynomial representation has some
negative coefficients.

These are distinct restrictions: exact reference monotonicity is proved
here for \(m^2\le8/3\), and residual monotonicity and reference domination
for \(m^2\le3\). All declared candidate endpoints must satisfy both when
all three facts are used.

## 3. Stronger error budget from the endpoint tail envelope

This is a new analytic budget lemma. It follows from the unsplit error
inequality and the complete endpoint certificate, not from an assumed
finite error norm. Define

\[
F(x)=-2x-4\log(1-x/2),\qquad F'(x)=x/\delta(x),
\]

so \(F(0)=0\) and \(F\ge0\) on \([0,2)\). Let
\(R_{\rm tot}=\int_0^\infty R\) and
\(Q_T=\int_0^T b h_e\).

The integrated endpoint transform gives the stronger time-dependent
envelope

\[
\Psi(b(t))\le M_{\rm true}e^{-(P-A(t))/4}
\le\Psi(B)e^{-(P-A(t))/4}.
\]

Define \(\beta(t)\in[0,B]\) by equality on the right. Then \(b\le\beta\),
\(\beta\) is increasing, and
\(\beta'=U^2\beta/(4\delta(\beta))\). Hence, using monotonicity of
\(x/\delta(x)\),

\[
J_T:=\int_0^T U^2\frac{b^2}{\delta(b)^2}
\le4\bigl[F(\beta(T))-F(\beta(0))\bigr]
\le4F(B).
\tag{2}
\]

Now multiply the original error inequality by \(b/\delta(b)\) and
integrate only to finite \(T\). The chain rule and Cauchy-Schwarz give

\[
F(b(T))-F(b(0))+Q_T
\le\sqrt{J_TQ_T}+\frac{B}{\delta(B)}R_{\rm tot}.
\]

Dropping the nonnegative terminal \(F(b(T))\), and using \(b(0)\le r\),
gives

\[
Q_T\le C+\sqrt{4F(B)Q_T},\qquad
C=F(r)+\frac{B}{\delta(B)}R_{\rm tot}.
\]

Solve this quadratic inequality before passing to infinite time:

\[
\sqrt Q\le\sqrt{F(B)}+\sqrt{F(B)+C},
\qquad Q=\int_0^\infty b h_e<\infty.
\tag{3}
\]

There is no circular assumption \(Q<\infty\). The finite-time estimates
hold first, then monotone convergence gives (3). Finally,
\(A_0(u-w_2)^2\le b h_e\), so

\[
\int_0^\infty A_0(u)^2
\le\left[\sqrt P+\sqrt{F(B)}+\sqrt{F(B)+C}\right]^2.
\tag{4}
\]

This avoids the severe \(\delta(B)^{-2}P\) loss of the elementary uniform
error-cap estimate. The term \(B R_{\rm tot}/\delta(B)\) and the logarithm
still diverge as \(B\uparrow2\); the proof requires an explicit cap
strictly below two.

## 4. Uniform Galerkin and NS-B implications

The second reference is supported in squared frequencies at most six.
For every cutoff \(N\ge\sqrt6\), it lies in the Galerkin range. Projection
can only decrease \(\chi(h)\) and the residual norm, and the same bilinear
inequalities apply to the projected output. The identical endpoint
envelope, \(\beta\), and finite-time estimates (2)--(4) therefore hold
uniformly for all these Galerkin solutions. The reference's known future
tail is not an assumption about the unknown solution's future.

If (4) is bounded by an explicit exponent \(K\), the existing
energy/enstrophy argument gives

\[
\sup_{N\ge\sqrt6}\int_0^\infty\|\nabla u_N\|_2^4
\le\frac{E_0X_0}{2}e^K,
\quad E_0=\|aW+h\|_2^2,\quad X_0=\|\nabla(aW+h)\|_2^2.
\]

For \(N<\sqrt6\), energy directly gives
\(\int\|\nabla u_N\|_2^4\le N^2E_0^2/2\le3E_0X_0\), using mean-zero
Poincare. Thus any \(K\) with \(e^K\ge6\) absorbs these small cutoffs
into the same final formula. No low-cutoff Fourier budget is asserted.

The estimates establish the usual uniform continuation and compactness
bounds without assuming an a priori uniform high Sobolev estimate.
The data-dependent product \(E_0X_0\) remains in the NS-B bound; a
critical-norm perturbation ball does not bound it uniformly.

## 5. Concrete enlarged regions

Every row below is independently certified. Here \(K\) is the integer
exponent in the uniform bound
\(\sup_N\int\|\nabla u_N\|_2^4\le(E_0X_0/2)e^K\).

| Amplitude interval | Critical perturbation radius | Error cap \(B\) | \(K\) |
|---|---|---|---|
| \(|a|\le1\) | \(\chi(h)\le1/10\) | \(9/5\) | 59 |
| \(|a|\le8/7\) | \(\chi(h)\le1/100\) | \(8/5\) | 49 |
| \(|a|\le11/10\) | \(\chi(h)\le1/40\) | \(6/5\) | 29 |
| \(|a|\le9/8\) | \(\chi(h)\le1/100\) | \(6/5\) | 30 |
| \(|a|\le1\) | \(\chi(h)\le2/25\) | \(6/5\) | 25 |

The last three rows give additional amplitude/radius tradeoffs or better
bounds on subsets. These classes form a union of separately specified
regions. Taking the amplitude maximum from one row and the perturbation
maximum from another would create an unsupported rectangle.

The two headline endpoint margins are positive: for \((m,r)=(1,1/10)\),
the certified endpoint upper bound is approximately \(0.7309181191561585\),
with margin above it to the barrier lower bound greater than
\(0.0009072683\); for \((m,r)=(8/7,1/100)\), they are approximately
\(0.7147619982450544\) and \(0.0041643443\). Exact rational comparisons,
not these displayed decimals, establish the strict inequalities.

The history-sensitive budget also improves the common exponent for the
two previously established classes from 81 to 24. This uses their
separately known endpoint cap \(B=1\), then bounds their inputs jointly
by \(r\le1/20\), \(P\le213/25\), \(R_{\rm tot}\le3/20\).
The resulting Fourier-budget upper bound is approximately
\(23.368970538848817\). It does not certify the unproved combined rectangle
\(|a|\le11/10,\chi(h)\le1/20\): the cap-one premise was established only
for the original union.

## 6. Exact-data audit and strict class comparison

The new `enlarged_endpoint_adversarial_static.py` returned
`PASS_INDEPENDENT_ENLARGED_ENDPOINT_AUDIT`. It independently rebuilt all
seven tested amplitude profiles from the nine-field residual formula,
recomputed their exact \(P_2\) and residual moments, and checked every
passing endpoint comparison using a separate 100-term positive Taylor
exponential enclosure. The residual series tail is valid because the
nonnegative time profile gives decreasing moments for all orders; beyond
term 25 their ratios are bounded by \(c/26<1\).

For all five new budgets and the old-class corollary, the checker
independently bounded \(F\) using a direct 300-term atanh series without
the budget program's dyadic normalization. It verified the saved outward
radical bounds by exact squaring inequalities and the final strict
integer-budget comparisons. The six exponents are
\(59,49,29,30,25\) for the rows above and 24 for the old union.

The audit is bound to these frozen sources:

- `enlarged_endpoint_gate_result.json`:
  `661a9c1c426b59e5f9d0a9477628f03f3b96a1ec8f86dc69f4c7d5e919ae4f0f`.
- `enlarged_budget_static_result.json`:
  `e809816968d76d780fbfa3a0ba145f73bcc317e2c79f8f930d243b358cb503c5`.

Complete check results are in
`enlarged_endpoint_adversarial_static_result.json`.

The root-owned class-inclusion source was also read and independently
checked mathematically, without executing it. Nonuniqueness of the
decomposition \(aW+h\) is accounted for. The old amplitude class is
contained in the new \((8/7,1/100)\) class and the old radius class in the
new \((1,1/10)\) class. Two explicit witnesses prove strict containment
relative to both named old classes:

1. For
   \(u_0=(\sin z,\sin x+(1/5)\sin2x,\sin y)\), disjoint Fourier supports
   give \(\chi(u_0-bW)=3|1-b|+1/10\). Thus no allowed value of \(b\)
   puts this field in either old radius, whereas the new unit-amplitude
   class includes it.
2. For \(u_0=(8/7)W\), the minimum critical distances to the old amplitude
   segments \(|b|\le11/10\) and \(|b|\le1\) are respectively
   \(9/70>1/100\) and \(3/7>1/20\). The field belongs to the new
   amplitude class.

Both witnesses are smooth, real, mean-zero and divergence-free. These
comparisons concern the two explicit old analytic classes, not every
possible earlier sufficient criterion or a novelty claim.

The enlarged endpoint gate establishes regularity and the stated uniform
NS-B bounds. It does not extend the separate Picard-convergence theorem
to these larger regions. Failed upper gates do not establish blow-up or
iteration divergence. No unresolved flaw was found within the stated
five regions; unrestricted-data regularity remains outside this result.
