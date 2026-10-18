# A concrete convergent residual-correction sequence for the unit W datum

This is a private local certificate for viscosity one on the normalized \(2\pi\)-periodic three-dimensional torus, with the single initial datum
\[
W=(\sin z,\sin x,\sin y).
\]
There is no initial perturbation in this certificate.

**Result.** The ordinary global Picard sequence
\[
w_0=e^{t\Delta}W,\qquad
w_{n+1}=e^{t\Delta}W-\mathcal Q(w_n,w_n)
\]
converges in \(\mathcal X=L^2(0,\infty;A_0)\) to its global smooth Navier–Stokes solution. In particular, the repeated residual heat correction is proved to converge here; it is not merely a decreasing-defect observation. The explicit certified rate is
\[
\boxed{\|w_n-u\|_{\mathcal X}\le14(3/4)^{n-3}\quad(n\ge3).}       \tag{1}
\]
A block invariant set around \(w_3\) also proves
\[
\boxed{\int_0^\infty A_0(u)^2\,dt<8.}                          \tag{2}
\]
The same construction for every solution Galerkin cutoff containing \(w_3\), followed by a separate energy argument for smaller cutoffs, gives
\[
\boxed{\sup_M\int_0^\infty
       \|\nabla u_M(t)\|_2^4\,dt\le\frac98e^8.}                 \tag{3}
\]
The claim is for this datum. No perturbation radius, arbitrary-data theorem, or universal success of residual correction is inferred.

## 1. Reference, defects, and exact iteration

Let \(B_s(a,b)=\tfrac12(B(a,b)+B(b,a))\), where \(B\) includes the Leray projection, and put
\[
\mathcal Q(a,b)(t)=\int_0^t e^{(t-s)\Delta}B_s(a,b)(s)\,ds,
\qquad
\mathcal H r(t)=\int_0^t e^{(t-s)\Delta}r(s)\,ds.
\]
Here \(A_0(f)=\sum_{k\ne0}|\widehat f_k|\), with Euclidean norm on each complex Fourier vector. The space \(\mathcal X\) consists of real, mean-zero, divergence-free fields with this time-integrated norm.
For a smooth reference with initial value \(W\),
\[
r_j=w_j'-\Delta w_j+B_s(w_j,w_j),\qquad
z_{j+1}=\mathcal H r_j,
\]
and hence
\[
w_j-z_{j+1}=e^{t\Delta}W-\mathcal Q(w_j,w_j)=w_{j+1}.            \tag{4}
\]
All finite stages retain the same initial value. The new residual is exactly
\[
r_{j+1}=-2B_s(w_j,z_{j+1})+B_s(z_{j+1},z_{j+1}).              \tag{5}
\]

The second reference, proved in the earlier second-correction notes, is
\[
w_2=(e^{-t}-L_1/2)W-te^{-2t}F+L_5C-L_6H,                     \tag{6}
\]
where \(F=B(W,W)\), \(H=B(F,F)\),
\[
C=\tfrac12(\sin z\cos2y,\sin x\cos2z,\sin y\cos2x),
\]
and
\[
\begin{split}
L_1&=[e^{-t}-(1+2t)e^{-3t}]/4,\\
L_5&=[e^{-5t}+(2t-1)e^{-3t}]/4,\\
L_6&=[(2t^2-2t+1)e^{-4t}-e^{-6t}]/4.
\end{split}
\]
Its disjoint Fourier shells give the exact amplitude
\[
\begin{split}
U_2=A_0(w_2)
 &=3e^{-t}+3te^{-2t}-\tfrac32(L_1-L_5)+\sqrt3L_6\\
 &=\tfrac{21}{8}e^{-t}+3te^{-2t}
       +\tfrac32te^{-3t}+\tfrac38e^{-5t}+\sqrt3L_6.             \tag{7}
\end{split}
\]
All \(L_j\) are nonnegative heat convolutions. Replacing \(\sqrt3\) by rational lower and upper enclosures in (7) therefore gives valid nonnegative amplitude bounds. The newly computed squared-norm enclosure is
\[
\|w_2\|_{\mathcal X}^2=6.569765510001933\ldots.                 \tag{8}
\]
This uses the actual \(U_2\), not its older upper profile \(U_1=3e^{-t}+3te^{-2t}\).

The center for the successful certificate is
\[
v=w_3=w_2-z_3,\qquad z_3=\mathcal H r_2.
\]
The full signed \(z_3\) has 342 modes and support \(|k|^2\le24\), so \(v\) also has support there. Its full next defect is
\[
z_4=\mathcal H[-2B_s(w_2,z_3)+B_s(z_3,z_3)].                  \tag{9}
\]
The fourth-defect calculation constructs the signed mixed term exactly and bounds the complete quadratic term. Its global enclosure is
\[
\|z_4\|_{\mathcal X}\le0.007982620403.
\]
More useful here are its block enclosures, including every earlier forcing time. Low output modes are retained in (9); there is no false lower frequency cutoff after correction.

## 2. A block heat kernel that retains earlier forcing

The incompressible symmetric Fourier estimate proved earlier implies
\[
\|P_k[(k\cdot a)b+(k\cdot b)a]\|\le |k|\,|a|\,|b|.
\]
Consequently
\[
\|\mathcal Q(a,b)\|_{\mathcal X}\le
c\|a\|_{\mathcal X}\|b\|_{\mathcal X},\qquad
c=\frac1{2\sqrt2}.                                           \tag{10}
\]
It holds for zero-extended time-block inputs and after Fourier output projections.

Partition time into finitely many intervals \(I_i=[t_{i-1},t_i]\), with the last endpoint infinity. For \(j<i\), write
\[
g=t_{i-1}-t_j\ge0,\qquad
\Delta=t_i-t_{i-1},
\]
and define
\[
K(g,\Delta)=\sup_{\lambda\ge1}
 \sqrt{e^{-2\lambda g}(1-e^{-2\lambda\Delta})},                 \tag{11}
\]
where \(1-e^{-2\lambda\infty}=1\). Choose \(k_{ij}\ge K(g,\Delta)\) and put \(k_{ii}=1\). Then
\[
\boxed{\|\mathcal Q(a,b)\|_{\mathcal X(I_i)}
 \le c\sum_{j\le i}k_{ij}
        \|a\|_{\mathcal X(I_j)}\|b\|_{\mathcal X(I_j)}.}       \tag{12}
\]

For completeness, split the bilinear forcing into its source intervals. For \(s\in I_j\), \(j<i\), a nonzero Fourier mode with \(\lambda=|k|^2\) has output-block kernel norm
\[
\left(\int_{I_i}|k|^2e^{-2|k|^2(t-s)}dt\right)^{1/2}
 =\frac1{\sqrt2}e^{-\lambda(t_{i-1}-s)}
          \sqrt{1-e^{-2\lambda\Delta}}
 \le\frac1{\sqrt2}K(g,\Delta).
\]
Minkowski in the causal source integral, followed by the symmetric Fourier bound and Cauchy–Schwarz in source time, proves each term of (12). The diagonal term uses the complete future kernel norm \(1/\sqrt2\). This argument sums all source blocks and does not restart the error at block boundaries. Only the universal torus gap \(|k|\ge1\) is used.

For finite \(g,\Delta>0\), the unrestricted maximum over \(\lambda\ge0\) is
\[
\max_{\lambda\ge0}e^{-2\lambda g}(1-e^{-2\lambda\Delta})
 =\frac{\Delta}{g+\Delta}
       \left(\frac g{g+\Delta}\right)^{g/\Delta}.              \tag{13}
\]
The maximizer is
\(\lambda_*=\log(1+\Delta/g)/(2\Delta)\).
If \(\lambda_*\le1\), the constrained maximum in (11) is instead attained at \(\lambda=1\). The checker uses outward rational bounds for the minimum of the unrestricted expression and \(e^{-2g}\), and uses the constrained value whenever \(\lambda_*\le1\) is itself certified. For \(g=0\) it uses one; for \(\Delta=\infty\) it uses \(e^{-g}\). Thus no floating point maximization enters the proof.

## 3. Componentwise invariant radii and convergence

Take certified nonnegative bounds
\[
V_i\ge\|v\|_{\mathcal X(I_i)},\qquad
d_i\ge\|\mathcal H r_v\|_{\mathcal X(I_i)}
       =\|z_4\|_{\mathcal X(I_i)}.
\]
Let \(\bar c\ge c\), and suppose positive radii \(\rho_i\) satisfy
\[
\boxed{d_i+\bar c\sum_{j\le i}
 k_{ij}(2V_j\rho_j+\rho_j^2)\le\rho_i.}                       \tag{14}
\]
Then \(T(f)=e^{t\Delta}W-\mathcal Q(f,f)\) preserves the closed product set
\[
\mathcal B=\{f:\|f-v\|_{\mathcal X(I_i)}\le\rho_i
                       \text{ for every }i\}.
\]
Indeed, \(T(v)-v=-z_4\), and (12) bounds the remaining linear and quadratic errors by (14).

For \(f,g\in\mathcal B\), define
\[
A_{ij}=2\bar c\,k_{ij}(V_j+\rho_j)\quad(j\le i),
\qquad A_{ij}=0\quad(j>i).                                   \tag{15}
\]
The vector of block distances obeys
\[
\big(\|T(f)-T(g)\|_{\mathcal X(I_i)}\big)_i
 \le A\big(\|f-g\|_{\mathcal X(I_i)}\big)_i.                   \tag{16}
\]
If each diagonal entry is below one, the finite nonnegative lower triangular matrix \(A\) has spectral radius below one. For any scalar \(q\) between that spectral radius and one, the positive vector
\(p=(qI-A)^{-1}\mathbf1\) satisfies \(Ap=qp-\mathbf1\le qp\).
The norm \(\max_i\|f\|_{\mathcal X(I_i)}/p_i\) is an equivalent complete norm on \(\mathcal X\), and (16) makes \(T\) a strict contraction on \(\mathcal B\).

This proves existence, uniqueness within the invariant set, and convergence of the actual global iteration. Starting from \(v=w_3\), the block error satisfies
\[
\big(\|T^n(v)-u\|_{\mathcal X(I_i)}\big)_i\le A^n\rho,
\quad
\boxed{\|T^n(v)-u\|_{\mathcal X}\le\|A^n\rho\|_{\ell^2}.}       \tag{17}
\]
These are complete time histories. There is no sequential solution restart, altered equation, or deletion of the low-mode residual.

## 4. The passing 17-block certificate

The partition endpoints are
\[
0,\tfrac1{16},\tfrac18,\tfrac3{16},\tfrac14,\tfrac5{16},
\tfrac38,\tfrac7{16},\tfrac12,\tfrac58,\tfrac34,\tfrac78,
1,\tfrac54,\tfrac32,2,3,\infty.
\]
The checker takes
\[
\bar c=\frac{707106781187}{2000000000000}>
\frac1{2\sqrt2},\qquad
V_i=\overline{\|w_2\|_{\mathcal X(I_i)}}+
       \overline{\|z_3\|_{\mathcal X(I_i)}}.                   \tag{18}
\]
The first term uses exact moments of (7). For the second, it sums the square roots of the exact modewise squared moments on \(I_i\), so signs and cancellations within each Fourier mode's complete time profile are retained before taking a norm. The \(d_i\) are the complete \(z_4\) block caps from the fourth-defect certificate.

The following table is a decimal orientation only. The exact rational \(V_i,d_i,\rho_i,k_{ij}\) and every row inequality are stored in *refined_W_block_gate_result.json* and are the operands of the proof.

| Block | Interval | \(V_i\) | \(d_i\) | \(\rho_i\) |
|---:|---|---:|---:|---:|
| 1 | 0–1/16 | 0.7486328552 | 0.000000516852 | 0.000001098198465 |
| 2 | 1/16–1/8 | 0.7410434914 | 0.000008978043 | 0.00002008291883 |
| 3 | 1/8–3/16 | 0.7274003727 | 0.000042046738 | 0.0001088543215 |
| 4 | 3/16–1/4 | 0.7090684209 | 0.000113437553 | 0.0003508934218 |
| 5 | 1/4–5/16 | 0.6871426313 | 0.000226498147 | 0.0008460415750 |
| 6 | 5/16–3/8 | 0.6625225410 | 0.000375386878 | 0.001694851418 |
| 7 | 3/8–7/16 | 0.6359596378 | 0.000549100632 | 0.002983786911 |
| 8 | 7/16–1/2 | 0.6080868880 | 0.000735155992 | 0.004771785083 |
| 9 | 1/2–5/8 | 0.7992626675 | 0.001425994865 | 0.01157548950 |
| 10 | 5/8–3/4 | 0.7176358846 | 0.001867485083 | 0.02185182746 |
| 11 | 3/4–7/8 | 0.6385494503 | 0.002188999518 | 0.03441503923 |
| 12 | 7/8–1 | 0.5641942910 | 0.002381351978 | 0.04789879375 |
| 13 | 1–5/4 | 0.6588725221 | 0.003399619248 | 0.08621836359 |
| 14 | 5/4–3/2 | 0.5018275018 | 0.003162118096 | 0.1170003465 |
| 15 | 3/2–2 | 0.4744112914 | 0.003461463471 | 0.1658877640 |
| 16 | 2–3 | 0.3083506673 | 0.002575157518 | 0.1811044031 |
| 17 | 3–infinity | 0.1072564237 | 0.000976407917 | 0.1138901612 |

Every row of (14) passes in exact rational arithmetic. The minimum quadratic discriminant is greater than \(0.18203\). The matrix satisfies
\[
\max_i A_{ii}=0.5733491592406857\ldots<3/5.                    \tag{19}
\]
No refined diagonal estimate is required for this successful certificate.

The exact radius vector gives
\[
\sum_i\rho_i^2=0.098536791493048\ldots,
\]
and, most usefully,
\[
\boxed{\sum_i(V_i+\rho_i)^2
 =7.508934101585295\ldots<8.}                                 \tag{20}
\]
This proves (2).

The independent finite rational rate calculation sets \(q=3/4\), solves
\[
p=(qI-A)^{-1}\mathbf1,\qquad
\alpha=\max_i(\rho_i/p_i),
\]
and verifies
\[
\alpha^2\sum_i p_i^2=183.535391\ldots<196.
\]
Since \(A^n\rho\le\alpha q^np\), (17) yields
\[
\|T^n(w_3)-u\|_{\mathcal X}\le14(3/4)^n\quad(n\ge0).
\]
The exact result is stored in *W_picard_rate_static_result.json*. Relation (4) translates it into (1). In particular, the standard heat-seed iteration converges after its finite prefix; the theorem does not require a different preconditioned map.

## 5. Every Galerkin cutoff and the NS-B bound

For \(M\ge\sqrt{24}\), \(P_Mv=v\). The projected map
\[
T_M(f)=e^{t\Delta}P_MW-P_M\mathcal Q(f,f)
\]
has center defect \(-P_Mz_4\). Projection contracts Fourier amplitude, so the same block bounds, invariant radii, and contraction matrix apply. This constructs the Galerkin solution inside the invariant set without assuming the desired full solution beforehand. In particular,
\[
\int_0^\infty A_0(u_M)^2dt<8\qquad(M\ge\sqrt{24}).             \tag{21}
\]

Write \(X_M=\|\nabla u_M\|_2^2\), \(Y_M=\|\Delta u_M\|_2^2\). The established coefficient-one \(H^1\) estimate gives
\[
X_M'+Y_M\le A_0(u_M)^2X_M,
\]
and therefore
\[
\sup_tX_M(t)\le X_0e^8.
\]
The energy identity gives \(\int_0^\infty X_M\le E_0/2\). Here \(E_0=X_0=3/2\), so
\[
\int_0^\infty X_M^2
 \le \frac{E_0X_0}{2}e^8=\frac98e^8
 \qquad(M\ge\sqrt{24}).                                     \tag{22}
\]

For smaller cutoffs the reference does not lie in the Galerkin space; (21) is not claimed. Instead, energy alone gives
\[
\sup_tX_M\le M^2E_0\le24(3/2)=36,\qquad
\int_0^\infty X_M\le3/4.
\]
Thus \(\int X_M^2\le27\). Since \(e^8>24\), this is at most the bound in (22), proving (3) for all cutoffs.

Uniform \(H^1\) and integrated \(H^2\) estimates give the usual global strong limit and smooth continuation for the smooth datum \(W\). The resulting solution is the mild fixed point above; local uniqueness and continuation identify the two constructions. The new analytic exponent eight applies to this exact datum. It does not replace earlier estimates on different amplitude or perturbation regions, nor a sharper earlier bound obtained under a different already-certified scope.

## 6. Bounded failed tests around w2, and their limitation

Before using \(w_3\), the checker examined three declared partitions with the exact \(U_2\) and the same heat-memory kernel:

- The 17-block partition used above.
- A 43-block partition: steps \(1/16\) through time 2, then \(1/8\) through time 3, then endpoints 4, 6, infinity.
- A 17-block partition: steps \(1/8\) through time 1, then \(1/4\) through time 2, then endpoints \(5/2\), 3, 4, 6, infinity.

For each, the center is \(w_2\) and the complete defect is \(z_3\). A second version replaces the diagonal linear coefficient \(2\bar cV_i\) by its minimum with the independently justified short-interval bound
\[
(17/25)\,U_1(t_{i-1})\sqrt{|I_i|};
\]
the diagonal quadratic coefficient remains \(\bar c\). This optional refinement is not used in Sections 3–5.

| Declared partition | First failed interval | Discriminant upper, original diagonal | Discriminant upper, refined diagonal |
|---|---|---:|---:|
| Original 17 | [3/4,7/8] | −0.0429565620 | −0.0249237195 |
| Fine 43 | [7/8,15/16] | −0.0714404246 | −0.0152606224 |
| Alternate 17 | [5/8,3/4] | −0.0354475921 | −0.0228223759 |

These failures refer to the displayed positive majorants. The checker propagates a lower bound on the smallest admissible preceding radii; a negative upper discriminant then precludes radii for that declared system. It does not prove divergence, failure of all partitions, or a lower bound on the actual defect equal to a conservative cap. In fact the next exact iterate \(w_3\) passes, so the original sequence starting from \(w_2\) converges as well.

## 7. Static verification, source binding, and scope

*refined_W_block_gate_check.py* is a new self-contained standard-library checker. It reads earlier JSON only as data, rebuilds the 342 signed modewise squares, and verifies every complete squared moment against its saved exact rational value. Those modes reduce to 19 distinct squared-norm time profiles.

All integrals use
\[
\int_t^\infty s^m e^{-\lambda s}ds
=e^{-\lambda t}\sum_{j=0}^m
 \frac{m!}{j!\lambda^{m-j+1}}t^j.
\]
The endpoint exponential enclosures use alternating rational Taylor bounds on \([0,1]\), followed by outward-rounded repeated squaring. Square roots use rational integer-square-root enclosures. Kernel logarithms use a finite positive series with an explicit geometric remainder. The sequential quadratic radii and all final inequalities use rational arithmetic. Decimal output is explanatory only.

Final bound sources are frozen with these SHA-256 values:

| Artifact | SHA-256 |
|---|---|
| third_heat_defect_static_result.json | b84689c2bd8cb746fd37a57cac76d83d7ca4f78ee401fb5a68621584c41003d0 |
| fourth_heat_defect_static_result.json | ab69f2fc5177e939cc985f5cc67f82d91990568a6338ee934c376840317465ec |
| refined_W_block_gate_result.json | 17b578bd37cb971cfbd21a58749bc95a94d62e126fcced0fe21bb3b85a3420b0 |

The block checker reports **PASS_REFINED_CONVERGENCE_CERTIFICATE**. The separate rate checker reports **PASS_UNIT_W_PICARD_EXPLICIT_GEOMETRIC_RATE**. An independent data-only adversarial checker verifies the row inequalities, heat-memory factors, matrix formulas, source binding, and budget.

There was no numerical evolution, quadrature, old-script execution or import, external search, upload, or change to earlier-stage artifacts. This is an analytic proof supported by finite exact arithmetic, not verification inside a formal proof assistant. The successful scope is the unit \(W\) datum with zero perturbation; extension to new data needs a new explicitly passing certificate.
