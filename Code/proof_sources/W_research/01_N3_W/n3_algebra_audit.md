# Independent algebra audit for the spherical N=3 reference

The script `n3_algebra_audit.py` independently reconstructs the Fourier
interactions in `n3_exact_algebra_intermediate.json`. It uses exact rational
arithmetic and performs no time integration. Run it with
`python -B n3_algebra_audit.py`; its compact certificate is
`n3_algebra_audit_result.json`. The certificate records the input SHA-256.

## Fourier and normalization conventions

Use normalized measure `(2*pi)^(-3) dx`. The stored real rational vector
`a_k` represents the imaginary Fourier coefficient `u_k = i*a_k`.
Consequently

\[
\frac{B(u,v)_k}{i}
=-\left(I-\frac{kk^\top}{|k|^2}\right)
\sum_{p+q=k}(a_p\cdot q)b_q.
\]

Parseval's squared norm is the sum of squared vector coefficients over all
signed modes. Conjugacy is `a_(-k) = -a_k`, and incompressibility is
`k dot a_k = 0`. The cutoff retains `0 < |k|^2 <= 9`.

For a diagonal basis pair, the audit computes `B(b_i,b_i)`. For `i < j`, it
computes `B(b_i,b_j)+B(b_j,b_i)`, without another factor two. The JSON
quadratic coefficients represent **positive** `P_3 B(v,v)`, so the ODE is

\[
x_\ell'=-s_\ell x_\ell-\sum_{i\le j}Q_{\ell ij}x_i x_j.
\]

## Exact checks

All 55 unordered basis-pair interactions are reconstructed independently.
The retained part is expanded in the supplied orthogonal basis using exact
inner products divided by its squared norms. The reconstructed field and
all 133 nonzero supplied ODE coefficients agree exactly. All closure
remainders vanish.

The ten basis vectors have squared-shell multiplicities
`1:1, 2:1, 5:2, 6:2, 8:1, 9:3`, with 108 signed retained modes. Their
normalized squared norms, in supplied order, are

\[
\left(\frac32,\frac34,\frac3{16},15,\frac18,12,6,1566,
\frac{87}{8},\frac32\right).
\]

Every basis vector has its declared shell, is mean-zero and divergence-free,
and satisfies conjugacy. All off-diagonal basis inner products vanish.
The labels W, C, D and F are also checked against independently reconstructed
definitions, including `C=B(W,W)` and `F=B(C,C)`.

The full exterior part is compared coefficient by coefficient for every
pair, including pairs whose exterior part is zero. There are 52 nonzero
exterior pairs and 648 signed exterior modes on squared shells

`10, 11, 13, 14, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30, 32, 33, 34`.

No exterior coefficients are discarded after the product is formed. The
audit finds no mismatches, duplicate keys, nonzero zero modes, divergence
errors or conjugacy errors. In particular,

\[
\|Q_3(B(W,D)+B(D,W))\|_2^2=\frac3{64}.
\]

Writing `m_l = ||b_l||_2^2`, the exact cubic polynomial
`sum_l m_l x_l sum_(i<=j) Q_(l i j) x_i x_j` vanishes identically. Thus the
reference obeys

\[
\frac{d}{dt}\sum_\ell m_\ell x_\ell^2
=-2\sum_\ell s_\ell m_\ell x_\ell^2.
\]

## Exact gradient norm of the first two basis vectors

For `W=(sin z,sin x,sin y)` and
`C=(sin y cos z,sin z cos x,sin x cos y)`, direct differentiation gives

\[
|\nabla(aW+bC)|_F^2
=\sum_{\mathrm{cyc}}(a\cos z-b\sin y\sin z)^2
+b^2\sum_{\mathrm{cyc}}\cos^2y\cos^2z.
\]

Cauchy-Schwarz bounds each first term by `a^2+b^2 sin^2 y`. Put
`p=sin^2 x`, `q=sin^2 y`, `r=sin^2 z`. Since
`pq+qr+rp <= p+q+r` on `[0,1]^3`,

\[
|\nabla(aW+bC)|_F^2
\le3a^2+b^2[3-(p+q+r)+(pq+qr+rp)]
\le3(a^2+b^2).
\]

Equality holds at the spatial origin. Therefore, for all real a and b,

\[
\boxed{\|\nabla(aW+bC)\|_\infty=\sqrt{3(a^2+b^2)}}.
\]

This proves the maximum rather than presuming its location. The audit
certifies the algebraic reference and residual data; numerical validation
and the a posteriori gate are separate checks.
