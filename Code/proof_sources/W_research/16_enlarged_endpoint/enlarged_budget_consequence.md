# Endpoint history gives a sharper NS-B budget

Private local analytic continuation, 2026-09-30 (Asia/Kolkata).
Earlier research files are read only. No previous checker is executed
or imported, and no time integration, quadrature, web access, or upload
is used.

The independently passing second-reference gates now certify, among
other classes, \(|a|\le1\), \(A_{-1}(h)\le1/10\), and
\(|a|\le8/7\), \(A_{-1}(h)\le1/100\), for initial data
\(aW+h\). The all-cutoff NS-B exponents below are respectively
\(59\) and \(49\). The two earlier classes also retain their proof
with a common exponent improved from \(81\) to \(24\).

The analytic gain is that an endpoint reference gate contains useful
information about the remaining reference energy. Retaining that
information gives a much smaller continuation budget than replacing
the critical error by its largest value throughout time. The result
applies to the full residual and to arbitrary Fourier support in the
initial perturbation. It does not assert ordinary Picard convergence
for these enlarged classes.

## 1. Assumptions and the endpoint envelope

Use the established viscosity-one periodic setting and let \(v\) be a
smooth divergence-free reference with full residual
\(r=v'+Av+B(v,v)\). For a Galerkin solution whose cutoff contains the
support of \(v\), put

\[
 b(t)=A_{-1}(u-v),\quad h(t)=A_1(u-v),\quad
 U(t)=A_0(v),\quad R(t)\ge A_{-1}(r).
\]

Projection can only reduce the initial critical error and the residual
norm. Assume \(b(0)\le r_0<2\), and write

\[
 P=\int_0^\infty U^2dt<\infty,\qquad
 R_* =\int_0^\infty Rdt<\infty,
 \qquad \Psi(s)=se^{-s/2}.
\]

The established unsplit comparison inequality is

\[
 b'+(1-b/2)h\le U\sqrt{bh}+R.                     \tag{1}
\]

Suppose the complete endpoint gate has been certified with some
\(0<B<2\):

\[
 M:=e^{P/4}\Psi(r_0)
   +\int_0^\infty e^{\frac14\int_t^\infty U^2}R(t)dt
       <\Psi(B).                                  \tag{2}
\]

The same argument applies when \(M\) in (2) is replaced by any
verified upper bound. The first-contact endpoint proof gives
\(b(t)<B\). In fact its partial-time calculation also gives

\[
 \Psi(b(t))\le
 \Psi(B)\exp\!\left(-\frac14\int_t^\infty U^2ds\right).    \tag{3}
\]

For completeness, the endpoint transform satisfies
\((\Psi(b))'\le U^2\Psi(b)/4+R\) below the barrier. Integrating
only through time \(t\), multiplying the resulting expression by
\(\exp(\frac14\int_t^\infty U^2)\), and extending the positive
forcing integral to infinity proves (3). Thus its direction does not
require a lower bound on the reference norm or on its past energy.

Define \(\beta(t)\) as the inverse of \(\Psi\) on its increasing
branch \([0,2)\), applied to the right side of (3). Then

\[
 b(t)\le\beta(t)\le B,\qquad
 \beta'=\frac{U^2\beta}{4(1-\beta/2)},\qquad
 \lim_{t\to\infty}\beta(t)=B.                     \tag{4}
\]

This is an explicit analytic envelope. No differential equation is
being numerically evolved.

## 2. A logarithmic dissipation bound

Set

\[
 \delta_s=1-s/2,\qquad
 F(s)=-2s-4\log(1-s/2).
\]

For \(0\le s<2\), \(F(0)=0\), \(F\ge0\), and
\(F'(s)=s/\delta_s\). Multiplying (1) by \(b/\delta_b\) gives

\[
 \frac{d}{dt}F(b)+bh
 \le\frac b{\delta_b}U\sqrt{bh}
       +\frac b{\delta_b}R.                       \tag{5}
\]

Because \(s/\delta_s\) increases and (4) bounds \(b\),

\[
 J:=\int_0^\infty U^2\frac{b^2}{\delta_b^2}dt
 \le\int_0^\infty U^2\frac{\beta^2}{\delta_\beta^2}dt
 =4\{F(B)-F(\beta(0))\}\le4F(B).                 \tag{6}
\]

The equality uses (4), equivalently
\(U^2dt=4\delta_\beta\,d\beta/\beta\). It remains valid if
\(U\) vanishes on some time intervals, either by the chain rule or
by absolute continuity; no strict monotonicity of \(\beta\) is
needed.

Let \(Q_T=\int_0^T bh\), and set

\[
 C=F(r_0)+\frac B{\delta_B}R_*.
\]

Integrate (5), drop the nonnegative final \(F(b(T))\), and use
Cauchy--Schwarz together with (6). For every finite \(T\),

\[
 Q_T\le C+\sqrt{4F(B)Q_T}.
\]

Solving this scalar quadratic first, and only then letting
\(T\to\infty\), proves without circular finiteness assumptions

\[
 \boxed{\sqrt{Q}\le\sqrt{F(B)}+\sqrt{F(B)+C},
 \qquad Q=\int_0^\infty bh\,dt.}                 \tag{7}
\]

Since \(A_0(u-v)^2\le bh\), Minkowski now gives the reference
comparison budget

\[
 \boxed{\int_0^\infty A_0(u)^2dt
 \le G:=\left[\sqrt P+\sqrt{F(B)}
                   +\sqrt{F(B)+C}\right]^2.}      \tag{8}
\]

All terms on the right have computable analytic upper bounds. The
initial term can also be bounded by
\(F(r_0)\le r_0^2/(2\delta_{r_0})\), from its integral formula,
if a rational estimate without a logarithm is preferred.

## 3. Comparison with a constant-barrier estimate

Discarding (3) and only using \(b\le B\) in (1) gives

\[
 \delta_B H\le r_0+R_*+\sqrt{BPH},\qquad
 H=\int_0^\infty hdt.
\]

Its optimal scalar quadratic bound is

\[
 \sqrt H\le
 \frac{\sqrt{BP}+\sqrt{BP+4\delta_B(r_0+R_*)}}
      {2\delta_B}.
\]

This improves a fixed Young-inequality split, but its dominant cost
scales like \(\delta_B^{-2}\). In (6), the coefficient instead
satisfies \(4F(B)=-8B-16\log\delta_B\), which diverges only
logarithmically as \(B\uparrow2\). The forcing contribution still
contains \(B R_*/\delta_B\); no claim is made to remove every
endpoint singularity.

## 4. Uniform cutoff consequences

The preceding proof starts with each finite Galerkin system; it does
not assume global smoothness of the full PDE. For every cutoff
containing the reference support, (8) and the established enstrophy
estimate imply

\[
 \sup_t\|\nabla u_K(t)\|_2^2\le X_0e^G,
 \qquad
 \int_0^\infty\|\nabla u_K(t)\|_2^4dt
       \le\frac{E_0X_0}{2}e^G,
\]

where \(E_0=\|u_0\|_2^2\) and
\(X_0=\|\nabla u_0\|_2^2\). The ordinary energy equality provides
\(\int\|\nabla u_K\|_2^2\le E_0/2\).

For the second corrected reference, the support radius is \(\sqrt6\).
Every lower cutoff \(K<\sqrt6\) is handled directly:

\[
 \sup_t\|\nabla u_K\|_2^2\le K^2E_0\le6X_0,
\]
\[
 \int_0^\infty\|\nabla u_K\|_2^4dt
 \le\frac{K^2E_0^2}{2}\le3E_0X_0.                 \tag{9}
\]

Here the last comparison uses the mean-zero Poincare inequality
\(E_0\le X_0\). Thus any reported integer exponent \(L\ge G\)
with \(e^L\ge6\) yields the original all-cutoff bound

\[
 \boxed{\sup_K\int_0^\infty\|\nabla u_K(t)\|_2^4dt
       \le\frac{E_0X_0}{2}e^L.}                   \tag{10}
\]

This lower-cutoff argument does not require the full initial energy to
be small within a critical-norm neighborhood. Uniform continuation
and the established Galerkin compactness argument then give a global
strong solution, and smooth initial data give global smoothness.

## 5. Concrete enlarged classes and integer exponents

The independently checked `enlarged_endpoint_gate_result.json` gives
passing endpoint gates for the following full smooth data classes,
with \(W=(\sin z,\sin x,\sin y)\). Every row means
\(|a|\le A\) and \(A_{-1}(h)\le r_0\); it covers the entire
amplitude segment and arbitrary Fourier support of \(h\).

| Amplitude bound \(A\) | Critical radius \(r_0\) | Error cap \(B\) | Budget \(G\), for orientation | Certified integer \(L\) in (10) |
|---:|---:|---:|---:|---:|
| \(1\) | \(2/25\) | \(6/5\) | \(24.40770\ldots\) | \(25\) |
| \(1\) | \(1/10\) | \(9/5\) | \(58.60140\ldots\) | \(59\) |
| \(11/10\) | \(1/40\) | \(6/5\) | \(28.06654\ldots\) | \(29\) |
| \(9/8\) | \(1/100\) | \(6/5\) | \(29.06160\ldots\) | \(30\) |
| \(8/7\) | \(1/100\) | \(8/5\) | \(48.27248\ldots\) | \(49\) |

The exact rational upper bounds for \(P\), \(R_*\), and the
endpoint margin are supplied by that frozen gate certificate. The new
budget checker substitutes those bounds into (8), encloses every
logarithm and square root outward, and checks \(G<L\) exactly.
The displayed decimals are explanatory and do not determine a sign.

For comparison, using the same certified input upper bounds, the
checker produces constant-barrier upper bounds above \(2379\) for
the second row and above \(730\) for the last row. These are not lower
bounds on actual solution budgets. Retaining the endpoint history is
therefore material to the resulting constants; this improvement is
not a choice of a slightly better Young-inequality parameter.

Relative to the earlier second-reference result, the unit-amplitude
critical radius doubles from \(1/20\) to \(1/10\). At radius
\(1/100\), the amplitude endpoint increases from \(11/10\) to
\(8/7\). The additional rows record useful tradeoffs. None of these
sufficient classes is asserted to be maximal, and failure of another
tested gate would not imply a singularity.

## 6. The previous two classes have a common exponent 24

The earlier proof established \(B=1\) for both

\[
 |a|\le11/10,\quad A_{-1}(h)\le1/100,
 \qquad\text{or}\qquad
 |a|\le1,\quad A_{-1}(h)\le1/20.
\]

Both have \(P<213/25\), \(r_0\le1/20\), and
\(R_*<3/20\). Formula (8) therefore gives the common strict bound

\[
 G\le\left[\sqrt{213/25}+\sqrt{F(1)}
       +\sqrt{F(1)+F(1/20)+3/10}\right]^2
       <24.                                       \tag{11}
\]

The new exact checker verifies (11). Consequently their old exponent
\(81\) can be replaced by \(24\) without changing those data
classes or re-running any earlier research program. The common caps
apply to the stated union, and do not certify the larger rectangle
\(|a|\le11/10\), \(A_{-1}(h)\le1/20\).

## 7. Verification and limits

The new standalone `enlarged_budget_static.py` reads the frozen gate
JSON as data and verifies 178 exact assertions. Logarithms are bounded
by normalizing their rational arguments to \([1,2]\) with powers of
two, then using the positive atanh series with 40 terms and an explicit
geometric remainder. Its series variable has magnitude at most
\(1/3\). Square roots are enclosed by integer square roots with
outward rational rounding at scale \(10^{30}\). These are finite
algebraic certificates, not numerical quadrature.

The output `enlarged_budget_static_result.json` records
`PASS_RETAINED_ENDPOINT_HISTORY_NSB_BUDGETS`, all exact upper bounds,
and all five integer exponents. The endpoint source SHA-256 is
`661a9c1c426b59e5f9d0a9477628f03f3b96a1ec8f86dc69f4c7d5e919ae4f0f`.
The analytic envelope and multiplier argument were also independently
checked by another team agent.

These conclusions enlarge an explicit regularity and NS-B data class.
They do not prove arbitrary-data three-dimensional regularity or
convergence of ordinary Picard iteration throughout the enlarged
classes. The separately certified ordinary-Picard neighborhood keeps
its own smaller radius and independent proof.
