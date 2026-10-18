# Fourth heat defect: exact mixed calculation and quadratic enclosure

Private local result at the fixed datum W=(sin z,sin x,sin y), viscosity one, normalized 2π-periodic torus. The initial perturbation is zero. All time computations below are exact exponential-polynomial algebra or rational enclosures of analytic expressions; no evolution equation or time quadrature was run.

## Result

For the previously constructed w2 and z3=H r2, let

\[
w_3=w_2-z_3,\qquad
r_3=-C(w_2,z_3)+B(z_3,z_3),\qquad z_4=\mathcal H r_3,
\]

where C(u,v)=B(u,v)+B(v,u) and H is the causal heat convolution. This is the ordinary Picard correction: w3=T_W(w2), and w4=w3-z4=T_W(w3). No Oseen or other modified iteration is substituted.

The new rigorous enclosure is

\[
\boxed{\|z_4\|_X\le\frac{7982620403}{10^{12}}
=0.007982620403},\qquad X=L^2_t A_0.
\]

The bound separates two different calculations:

| Contribution | Method | X-norm upper bound |
|---|---|---:|
| -H C(w2,z3) | Exact Fourier and time-polynomial construction; exact modewise squared moments | 0.007361535372 |
| H B(z3,z3) | Proven bilinear estimate applied to a rigorous pointwise bound for A0(z3) | 0.000621085031 |

The mixed field has **1,290 modes**, maximum squared wave number **50**, and **43,674 vector temporal terms**. The quadratic field was not expanded. In particular, the mode count describes the exact mixed contribution, not the whole z4 field.

## Exact mixed term and source provenance

`fourth_heat_defect_static.py` reconstructs W,F,G5,H and the four w2 time profiles from scratch. It reads only `third_heat_defect_static_result.json` as historical data, whose SHA-256 is

`b84689c2bd8cb746fd37a57cac76d83d7ca4f78ee401fb5a68621584c41003d0`.

The code forms -C(w2,z3), retaining signs and combining all contributions to each Fourier mode before applying its heat convolution. Every mixed mode is checked for zero initial value, divergence-free structure, and its exact forced heat identity. For each mode k, its integral squared norm is computed by

\[
\int_0^\infty t^m e^{-\rho t}\,dt=\frac{m!}{\rho^{m+1}}.
\]

Minkowski then gives the stated mixed X bound. The full signed mixed histories are stored in `fourth_heat_defect_static_result.json` under `mixed_mode_profiles`.

## Quadratic bound and retained time history

Let Z(t) be the previously certified pointwise upper profile for A0(z3(t)). Its complete profile has been combined before integrating its square. The exact rational moment has decimal value

\[
\int_0^\infty Z(t)^2\,dt
=0.0017566937482090424\ldots.
\]

The full rational value is in the JSON. It also improves the previous z3 upper enclosure to 0.041912930561. With c=1/(2√2),

\[
\|\mathcal H B(z_3,z_3)\|_X
\le c\|z_3\|_X^2
\le c\int Z^2.
\]

For a block [a,b], including b=∞, the causal kernel gives the stronger block bound

\[
\|\mathcal H B(z_3,z_3)\|_{X(a,b)}
\le c\left[
\int_0^a e^{-(a-s)}Z(s)^2\,ds+
\int_a^b Z(s)^2\,ds\right].
\]

The first term explicitly preserves all earlier forcing. The kernel estimate uses nonzero torus wave numbers |k|≥1. Taking the minimum with the global quadratic bound remains valid.

The JSON provides full z4 norm upper bounds on the same 17 blocks with endpoints

0, 1/16, 1/8, 3/16, 1/4, 5/16, 3/8, 7/16, 1/2, 5/8, 3/4, 7/8, 1, 5/4, 3/2, 2, 3, ∞.

Mixed block integrals use exact antiderivatives and outward rational bounds for endpoint exponentials. The quadratic blocks use the displayed causal history formula. Their rational caps are under `block_defect_bounds`; no time interval resets the defect to zero.

## A genuine single-step decrease, with its precise limitation

Comparing only two upper bounds would not establish that the actual defect decreased. Here an independent lower bound is available from the exact z3 mode moments:

\[
\sum_k\int_0^\infty |\widehat z_{3,k}(t)|^2\,dt
=\frac{2628329906545220056349884107745448088259779767}
{30074849493343772457434740248769174123315200000000}.
\]

Because the Fourier l1 norm dominates l2,

\[
\|z_3\|_X\ge\frac{9348419809}{10^{12}}
=0.009348419809
>0.007982620403\ge\|z_4\|_X.
\]

Thus **this single actual Picard defect decreases**. A valid ratio bound is

\[
\frac{\|z_4\|_X}{\|z_3\|_X}
\le\frac{7982620403}{9348419809}
<0.854.
\]

This is a result about one transition. It does not prove a contraction for all subsequent steps, convergence of the whole iteration, or a result for arbitrary initial data. Any initial perturbation must be retained separately in a subsequent comparison or contraction gate.

The standalone checker reports **PASS_FOURTH_HEAT_DEFECT_STATIC_CHECKS**. Its checks are finite exact arithmetic supporting the analytic estimates, not formal proof-system verification. All previous files remain unchanged.
