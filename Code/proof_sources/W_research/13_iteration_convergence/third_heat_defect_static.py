"""Fresh exact heat-defect calculation at a=1; no numerical evolution.

Fourier coefficients are i times rational-vector exponential polynomials.
All time integrals use exact factorial identities; no quadrature is used.
No historical source, scripts, modules, or result files are loaded.
"""

from fractions import Fraction as Q
from math import factorial, isqrt
from pathlib import Path
from functools import lru_cache
import json


ZERO = (Q(0), Q(0), Q(0))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def va(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vs(c, a):
    return tuple(c * x for x in a)


def norm2(a):
    return dot(a, a)


def qs(q):
    q = Q(q)
    return f"{q.numerator}/{q.denominator}"


def fa(u, v):
    out = dict(u)
    for k, x in v.items():
        out[k] = va(out.get(k, ZERO), x)
    return {k: x for k, x in out.items() if x != ZERO}


def fs(c, u):
    return {k: vs(c, x) for k, x in u.items() if c}


def B(u, v):
    raw = {}
    for p, x in u.items():
        for q, y in v.items():
            k = va(p, q)
            raw[k] = va(raw.get(k, ZERO), vs(-dot(x, q), y))
    out = {}
    for k, x in raw.items():
        lam = norm2(k)
        if not lam:
            assert x == ZERO
            continue
        val = va(x, vs(-dot(k, x) / lam, k))
        if val != ZERO:
            out[k] = val
    return out


def C(u, v):
    return fa(B(u, v), B(v, u))


def pa(*profiles):
    out = {}
    for profile in profiles:
        for key, value in profile.items():
            out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}


def ps(c, profile):
    return {key: c * value for key, value in profile.items() if c * value}


def pm(a, b):
    out = {}
    for (r, n), x in a.items():
        for (s, m), y in b.items():
            key = (r + s, n + m)
            out[key] = out.get(key, Q(0)) + x * y
    return {key: value for key, value in out.items() if value}


def p_integral(profile):
    return sum((value * Q(factorial(n), r ** (n + 1))
                for (r, n), value in profile.items()), Q(0))


def p_heat(profile, lam):
    lam = int(lam)
    out = {}
    for (rate, power), value in profile.items():
        if rate == lam:
            out = pa(out, {(rate, power + 1): value / (power + 1)})
            continue
        delta = lam - rate
        particular = {
            (rate, j): value * Q((-1) ** (power - j) * factorial(power),
                                 factorial(j) * delta ** (power - j + 1))
            for j in range(power + 1)
        }
        out = pa(out, particular, {(lam, 0): -particular[(rate, 0)]})
    return out


def p_operator(profile, lam):
    out = {}
    for (rate, power), value in profile.items():
        out = pa(out, {(rate, power): (lam - rate) * value})
        if power:
            out = pa(out, {(rate, power - 1): power * value})
    return out


def p_initial(profile):
    return sum((value for (rate, power), value in profile.items() if power == 0), Q(0))


def p_json(profile):
    return [{"rate": r, "power": n, "coefficient": qs(value)}
            for (r, n), value in sorted(profile.items())]


def sqrt_bounds(q, denominator=10**12):
    q = Q(q)
    assert q >= 0
    # Compute floor(sqrt(q)*denominator) using integer arithmetic only.
    n = isqrt((q.numerator * denominator**2) // q.denominator)
    lower = Q(n, denominator)
    upper = lower if lower**2 == q else Q(n + 1, denominator)
    assert lower**2 <= q <= upper**2
    return lower, upper


@lru_cache(maxsize=None)
def negative_exp_bounds(x):
    """Exact Taylor/geometric bounds for exp(-x), rounded outward rationally."""
    x = Q(x)
    if not x:
        return Q(1), Q(1)
    assert x > 0
    final_index = max(80, (3 * x.numerator + x.denominator - 1) // x.denominator + 40)
    partial = term = Q(1)
    for n in range(1, final_index + 1):
        term *= x / n
        partial += term
    omitted = term * x / (final_index + 1)
    ratio = x / (final_index + 2)
    assert 0 < ratio < 1
    positive_upper = partial + omitted / (1 - ratio)
    lo, hi = 1 / positive_upper, 1 / partial
    denominator = 10**90
    rounded_lo = Q((lo.numerator * denominator) // lo.denominator, denominator)
    rounded_hi = Q((hi.numerator * denominator + hi.denominator - 1) // hi.denominator, denominator)
    assert rounded_lo <= lo <= hi <= rounded_hi
    return rounded_lo, rounded_hi


def tail_polynomial(profile):
    """Return exact profile for integral_t^infinity profile(s) ds."""
    out = {}
    for (rate, power), value in profile.items():
        out = pa(out, {(rate, j): value * Q(factorial(power), factorial(j) * rate ** (power - j + 1))
                       for j in range(power + 1)})
    return out


def endpoint_bounds(profile, time):
    if time is None:
        return Q(0), Q(0)
    by_rate = {}
    for (rate, power), value in profile.items():
        by_rate[rate] = by_rate.get(rate, Q(0)) + value * time**power
    lo = hi = Q(0)
    for rate, value in by_rate.items():
        elo, ehi = negative_exp_bounds(rate * time)
        if value >= 0:
            lo += value * elo
            hi += value * ehi
        else:
            lo += value * ehi
            hi += value * elo
    return lo, hi


def main():
    checks = []

    def check(name, ok):
        assert ok, name
        checks.append(name)

    W = {}
    for axis, component in ((2, 0), (0, 1), (1, 2)):
        for sign in (-1, 1):
            k = tuple(sign if j == axis else 0 for j in range(3))
            W[k] = tuple(Q(-sign, 2) if j == component else Q(0) for j in range(3))
    F = B(W, W)
    G = C(W, F)
    G1 = {k: v for k, v in G.items() if norm2(k) == 1}
    G5 = {k: v for k, v in G.items() if norm2(k) == 5}
    H = B(F, F)
    check("G1_minus_W_over_two", G1 == fs(Q(-1, 2), W))
    check("G_split", G == fa(G1, G5))

    p = {(1, 0): Q(1)}
    q = {(2, 1): Q(1)}
    L1 = p_heat({(3, 1): Q(1)}, 1)
    L5 = p_heat({(3, 1): Q(1)}, 5)
    L6 = p_heat({(4, 2): Q(1)}, 6)
    d, e, f = ps(Q(1, 2), L1), L5, L6
    # r2 has the following signed field coefficients. Their absolute values
    # are nonnegative for a=1 since d<=p/8, and all p,q,d,e,f are nonnegative.
    terms = [
        ("F", F, pa(pm(d, d), ps(-2, pm(p, d))), -1),
        ("G", G, pm(q, d), 1),
        ("C_W_G5", C(W, G5), pm(pa(p, ps(-1, d)), e), 1),
        ("C_F_G5", C(F, G5), ps(-1, pm(q, e)), -1),
        ("C_W_H", C(W, H), ps(-1, pm(pa(p, ps(-1, d)), f)), -1),
        ("C_F_H", C(F, H), pm(q, f), 1),
        ("B_G5_G5", B(G5, G5), pm(e, e), 1),
        ("B_H_H", B(H, H), pm(f, f), 1),
        ("C_G5_H", C(G5, H), ps(-1, pm(e, f)), -1),
    ]
    residual = {}
    for name, field, profile, sign in terms:
        for k, vector in field.items():
            if k not in residual:
                residual[k] = [{}, {}, {}]
            for j in range(3):
                residual[k][j] = pa(residual[k][j], ps(vector[j], profile))
    residual = {k: vec for k, vec in residual.items() if any(vec)}
    defect = {k: [p_heat(component, norm2(k)) for component in vector]
              for k, vector in residual.items()}
    for k, vector in defect.items():
        for j in range(3):
            check("zero_initial_" + str(k) + "_" + str(j), p_initial(vector[j]) == 0)
            check("heat_identity_" + str(k) + "_" + str(j),
                  p_operator(vector[j], norm2(k)) == residual[k][j])
        check("solenoidal_" + str(k), pa(*(ps(k[j], vector[j]) for j in range(3))) == {})

    # Rigorous Minkowski: ||A0(z)||L2 <= sum_k ||z_k||L2.
    mode_records = []
    sum_lower = sum_upper = Q(0)
    groups = {}
    squared_profile_groups = {}
    for k, vector in sorted(defect.items()):
        squared_profile = pa(*(pm(component, component) for component in vector))
        integral_square = p_integral(squared_profile)
        profile_key = tuple(sorted(squared_profile.items()))
        squared_profile_groups[profile_key] = squared_profile_groups.get(profile_key, 0) + 1
        check("nonnegative_exact_square_integral_" + str(k), integral_square >= 0)
        lo, hi = sqrt_bounds(integral_square)
        sum_lower += lo
        sum_upper += hi
        groups[integral_square] = groups.get(integral_square, 0) + 1
        mode_records.append({"k": list(k), "lambda": qs(norm2(k)),
                             "exact_L2_time_squared": qs(integral_square),
                             "L2_time_norm_upper": qs(hi),
                             "imaginary_vector_profiles": [p_json(component) for component in vector]})

    # A pointwise-valid block-friendly upper profile retains positive causal
    # convolutions of each signed scalar coefficient instead of taking absolute
    # values of its canceling exponential expansion.
    pointwise_terms = []
    field_shell_upper = Q(0)
    for name, field, profile, sign in terms:
        absolute_profile = ps(sign, profile)
        shells = {}
        for k, value in field.items():
            lam = norm2(k)
            shells.setdefault(lam, []).append(norm2(value))
        for lam, squares in sorted(shells.items()):
            amplitude_upper = sum((sqrt_bounds(s)[1] for s in squares), Q(0))
            heat_profile = p_heat(absolute_profile, lam)
            square_integral = p_integral(pm(heat_profile, heat_profile))
            check("nonnegative_shell_profile_square_" + name + "_" + str(lam), square_integral >= 0)
            field_shell_upper += amplitude_upper * sqrt_bounds(square_integral)[1]
            pointwise_terms.append({
                "field": name, "lambda": qs(lam), "amplitude_rational_upper": qs(amplitude_upper),
                "absolute_residual_coefficient": p_json(absolute_profile),
                "causal_heat_profile": p_json(heat_profile),
                "exact_scalar_L2_squared": qs(square_integral),
                "positivity": "Scalar forcing is nonnegative by the signed factorization; causal heat convolution is nonnegative.",
            })

    # Exact low-mode history, not a fresh zero-error restart.
    J_source = pa(ps(Q(1, 2), pm(q, d)), ps(Q(1, 24), pm(e, f)))
    J = p_heat(J_source, 1)
    J5 = p_heat(ps(Q(1, 2), pm(q, d)), 1)
    J7 = p_heat(ps(Q(1, 24), pm(e, f)), 1)
    J5_bracket = {(0, 0): Q(1), (2, 0): Q(-1), (2, 1): Q(-2)}
    J5_closed = ps(Q(1, 128), pm({(1, 0): Q(1)}, pm(J5_bracket, J5_bracket)))
    check("fundamental_history_J5_plus_J7", J == pa(J5, J7))
    check("J5_closed_positive_square", J5 == J5_closed)
    third_reference_lower_P = Q(3)**2 * Q(7, 8)**2 / 2
    check("third_reference_global_ball_obstruction", third_reference_lower_P == Q(441, 128) > 2)
    for k, value in W.items():
        for j in range(3):
            check("fundamental_profile_" + str(k) + "_" + str(j),
                  defect[k][j] == ps(-value[j], J))
    J_squared = p_integral(pm(J, J))
    J_lo, J_hi = sqrt_bounds(J_squared)
    lower_defect_norm = 3 * J_lo

    # One requested finite grid only. Block norms preserve the complete z3
    # history; they do not reset its left endpoint value to zero.
    grid = [Q(0), Q(1, 16), Q(1, 8), Q(3, 16), Q(1, 4), Q(5, 16),
            Q(3, 8), Q(7, 16), Q(1, 2), Q(5, 8), Q(3, 4), Q(7, 8),
            Q(1), Q(5, 4), Q(3, 2), Q(2), Q(3), None]
    block_bounds = [{"left": qs(left), "right": "infinity" if right is None else qs(right),
                     "upper": Q(0), "sum_lower": Q(0)} for left, right in zip(grid, grid[1:])]
    for key, count in squared_profile_groups.items():
        antiderivative = tail_polynomial(dict(key))
        endpoints = [endpoint_bounds(antiderivative, t) for t in grid]
        for i, block in enumerate(block_bounds):
            integral_lo = endpoints[i][0] - endpoints[i + 1][1]
            integral_hi = endpoints[i][1] - endpoints[i + 1][0]
            assert integral_hi >= 0
            lo, hi = sqrt_bounds(max(Q(0), integral_lo))[0], sqrt_bounds(integral_hi)[1]
            block["upper"] += count * hi
            block["sum_lower"] += count * lo
    fundamental_tail = tail_polynomial(pm(J, J))
    fundamental_endpoints = [endpoint_bounds(fundamental_tail, t) for t in grid]
    for i, block in enumerate(block_bounds):
        fundamental_integral_lower = max(Q(0), fundamental_endpoints[i][0] - fundamental_endpoints[i + 1][1])
        block["defect_norm_lower_bound"] = qs(3 * sqrt_bounds(fundamental_integral_lower)[0])
        block["defect_norm_upper_bound"] = qs(block.pop("upper"))
        block["Minkowski_sum_lower_bound_not_defect_lower"] = qs(block.pop("sum_lower"))
        block["upper_decimal_for_orientation"] = float(Q(block["defect_norm_upper_bound"]))
    check("seventeen_requested_blocks", len(block_bounds) == 17)
    result = {
        "status": "PASS_THIRD_HEAT_DEFECT_STATIC_CHECKS",
        "scope": "Fresh self-contained Fourier construction and exact rational exponential-polynomial convolution/moments at a=1. No old code/data loaded, numerical time integration, quadrature, or unknown-flow simulation.",
        "definitions": {
            "z3": "z3=integral_0^t exp((t-s)Delta) r2(s) ds",
            "w3": "w3=w2-z3",
            "r2": "Exact residual of w2, with all nine spatial pair fields combined mode by mode before convolution.",
            "p_q_d_e_f": "p=e^-t,q=t e^-2t,d=L1/2,e=L5,f=L6",
        },
        "defect_support_count": len(defect),
        "max_squared_wavenumber": qs(max(norm2(k) for k in defect)),
        "global_X_norm_bound": {
            "formula": "||z3||_X <= sum_k sqrt(integral |z3_k(t)|^2 dt)",
            "rational_upper_bound": qs(sum_upper),
            "decimal_upper_for_orientation": float(sum_upper),
            "rational_lower_bound_on_Minkowski_sum_not_on_X": qs(sum_lower),
            "exact_radical_grouping": [{"count": count, "sqrt_rational": qs(value)} for value, count in sorted(groups.items())],
            "field_shell_triangle_upper_bound": qs(field_shell_upper),
        },
        "pointwise_A0_upper_profile": {
            "formula": "A0(z3(t)) <= sum of amplitude_rational_upper times causal_heat_profile over these terms.",
            "use": "Each heat profile is nonnegative despite signed coefficients in its exponential-polynomial expansion. This upper profile can be integrated on blocks using exact moments or rigorous exponential bounds.",
            "terms": pointwise_terms,
        },
        "fundamental_history": {
            "formula": "z3_fundamental=-J(t)W, J=H1(qd/2+ef/24)>=0",
            "J_source": p_json(J_source),
            "J_profile": p_json(J),
            "J_exact_L2_squared": qs(J_squared),
            "X_norm_lower_bound_from_fundamental": qs(lower_defect_norm),
            "decimal_lower_for_orientation": float(lower_defect_norm),
            "J5_closed_form": "exp(-t)[1-(1+2t)exp(-2t)]^2/128",
            "J7": "H1(ef/24)>=0",
        },
        "third_reference_global_ball_obstruction": {
            "w3_fundamental": "[exp(-t)-L1/2+J5+J7]W >= (7/8)exp(-t)W in scalar coefficient",
            "P3_lower_bound": qs(third_reference_lower_P),
            "consequence": "P3>2, so the sufficient simple global contraction-ball hypothesis 2c||w3||_X<1 with c=1/(2sqrt2) fails. This is not a failure of other criteria or a nonexistence result.",
        },
        "block_defect_bounds": block_bounds,
        "block_bound_method": "Group identical exact squared Fourier coefficient profiles, integrate each exactly via exponential-polynomial antiderivatives, enclose endpoint exponentials by rational positive Taylor/geometric tails, then apply modewise L2 Minkowski. Full causal history retained.",
        "all_defect_modes": mode_records,
        "number_of_assertions": len(checks),
        "assertions": checks,
    }
    destination = Path(__file__).with_name("third_heat_defect_static_result.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "assertions": len(checks),
                      "mode_count": len(defect), "max_squared_wavenumber": result["max_squared_wavenumber"],
                      "X_upper": qs(sum_upper), "X_upper_decimal": float(sum_upper),
                      "X_lower_fundamental": qs(lower_defect_norm), "X_lower_decimal": float(lower_defect_norm),
                      "field_shell_X_upper_decimal": float(field_shell_upper),
                      "pointwise_profile_term_count": len(pointwise_terms),
                      "block_upper_bounds": [b["upper_decimal_for_orientation"] for b in block_bounds],
                      "result_file": str(destination)}, indent=2))


if __name__ == "__main__":
    main()
