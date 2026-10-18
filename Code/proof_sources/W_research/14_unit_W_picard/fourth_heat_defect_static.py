"""Fresh exact mixed fourth heat defect and rigorous quadratic envelopes.

The previous z3 JSON is read-only data. No historical scripts are imported or
executed. All temporal calculations are rational exponential-polynomial
algebra and exact moments, with rational endpoint exponential enclosures.
"""

from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json


ZERO = (Q(0), Q(0), Q(0))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def va(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vs(c, v):
    return tuple(c * x for x in v)


def n2(v):
    return dot(v, v)


def qs(q):
    q = Q(q)
    return f"{q.numerator}/{q.denominator}"


def pa(*profiles):
    out = {}
    for p in profiles:
        for k, v in p.items():
            out[k] = out.get(k, Q(0)) + v
    return {k: v for k, v in out.items() if v}


def ps(c, profile):
    return {k: c * v for k, v in profile.items() if c * v}


def pm(a, b):
    out = {}
    for (r, n), x in a.items():
        for (s, m), y in b.items():
            key = (r + s, n + m)
            out[key] = out.get(key, Q(0)) + x * y
    return {k: v for k, v in out.items() if v}


def pi(profile):
    return sum((v * Q(factorial(n), r ** (n + 1)) for (r, n), v in profile.items()), Q(0))


def tail(profile):
    out = {}
    for (r, n), v in profile.items():
        for j in range(n + 1):
            key = (r, j)
            out[key] = out.get(key, Q(0)) + v * Q(factorial(n), factorial(j) * r ** (n - j + 1))
    return {k: v for k, v in out.items() if v}


def pjson(p):
    return [{"rate": r, "power": n, "coefficient": qs(v)} for (r, n), v in sorted(p.items())]


def parse_profile(entries):
    return {(x["rate"], x["power"]): Q(x["coefficient"]) for x in entries}


def B(u, v):
    raw = {}
    for p, x in u.items():
        for q, y in v.items():
            k = va(p, q)
            raw[k] = va(raw.get(k, ZERO), vs(-dot(x, q), y))
    out = {}
    for k, v in raw.items():
        lam = n2(k)
        if lam == 0:
            assert v == ZERO
            continue
        projected = va(v, vs(-dot(k, v) / lam, k))
        if projected != ZERO:
            out[k] = projected
    return out


def fa(u, v):
    out = dict(u)
    for k, x in v.items():
        out[k] = va(out.get(k, ZERO), x)
    return {k: v for k, v in out.items() if v != ZERO}


def vpa_into(target, key, value):
    if value != ZERO:
        target[key] = va(target.get(key, ZERO), value)


@lru_cache(maxsize=None)
def heat_monomial(rate, power, lam):
    if rate == lam:
        return (((rate, power + 1), Q(1, power + 1)),)
    delta = lam - rate
    out = {(rate, j): Q((-1) ** (power - j) * factorial(power),
                        factorial(j) * delta ** (power - j + 1)) for j in range(power + 1)}
    out[(lam, 0)] = -out[(rate, 0)]
    return tuple(out.items())


def vector_heat(profile, lam):
    out = {}
    for (rate, power), vec in profile.items():
        for key, coefficient in heat_monomial(rate, power, lam):
            vpa_into(out, key, vs(coefficient, vec))
    return {k: v for k, v in out.items() if v != ZERO}


def vector_operator(profile, lam):
    out = {}
    for (rate, power), vec in profile.items():
        vpa_into(out, (rate, power), vs(lam - rate, vec))
        if power:
            vpa_into(out, (rate, power - 1), vs(power, vec))
    return {k: v for k, v in out.items() if v != ZERO}


def vector_square(profile):
    entries = list(profile.items())
    out = {}
    for i, ((r, n), x) in enumerate(entries):
        for j in range(i, len(entries)):
            (s, m), y = entries[j]
            value = dot(x, y) * (1 if i == j else 2)
            if value:
                key = (r + s, n + m)
                out[key] = out.get(key, Q(0)) + value
    return {k: v for k, v in out.items() if v}


def sqrt_bounds(q, denominator=10**12):
    q = Q(q)
    assert q >= 0
    n = isqrt(q.numerator * denominator**2 // q.denominator)
    lo = Q(n, denominator)
    hi = lo if lo**2 == q else Q(n + 1, denominator)
    assert lo**2 <= q <= hi**2
    return lo, hi


@lru_cache(maxsize=None)
def expneg(x):
    x = Q(x)
    if not x:
        return Q(1), Q(1)
    assert x > 0
    last = max(80, (3 * x.numerator + x.denominator - 1) // x.denominator + 40)
    partial = term = Q(1)
    for n in range(1, last + 1):
        term *= x / n
        partial += term
    ratio = x / (last + 2)
    upper = partial + term * x / (last + 1) / (1 - ratio)
    lo, hi = 1 / upper, 1 / partial
    denominator = 10**90
    return (Q(lo.numerator * denominator // lo.denominator, denominator),
            Q((hi.numerator * denominator + hi.denominator - 1) // hi.denominator, denominator))


def endpoint(profile, t):
    if t is None:
        return Q(0), Q(0)
    rates = {}
    for (r, n), v in profile.items():
        rates[r] = rates.get(r, Q(0)) + v * t**n
    lo = hi = Q(0)
    for r, value in rates.items():
        a, b = expneg(r * t)
        lo += value * (a if value >= 0 else b)
        hi += value * (b if value >= 0 else a)
    return lo, hi


def round_upper(q, denominator=10**12):
    q = Q(q)
    return Q((q.numerator * denominator + q.denominator - 1) // q.denominator, denominator)


def main():
    source = Path(__file__).with_name("third_heat_defect_static_result.json")
    source_bytes = source.read_bytes()
    previous = json.loads(source_bytes)
    assert previous["status"] == "PASS_THIRD_HEAT_DEFECT_STATIC_CHECKS"
    z3 = {}
    for mode in previous["all_defect_modes"]:
        profile = {}
        for j, scalar_profile in enumerate(mode["imaginary_vector_profiles"]):
            for key, value in parse_profile(scalar_profile).items():
                v = tuple(value if i == j else Q(0) for i in range(3))
                vpa_into(profile, key, v)
        z3[tuple(mode["k"])] = profile

    W = {}
    for axis, component in ((2, 0), (0, 1), (1, 2)):
        for sign in (-1, 1):
            k = tuple(sign if i == axis else 0 for i in range(3))
            W[k] = tuple(Q(-sign, 2) if i == component else Q(0) for i in range(3))
    F = B(W, W)
    G = fa(B(W, F), B(F, W))
    G5 = {k: v for k, v in G.items() if n2(k) == 5}
    H = B(F, F)
    L1 = {(1, 0): Q(1, 4), (3, 0): Q(-1, 4), (3, 1): Q(-1, 2)}
    L5 = {(5, 0): Q(1, 4), (3, 0): Q(-1, 4), (3, 1): Q(1, 2)}
    L6 = {(4, 0): Q(1, 4), (4, 1): Q(-1, 2), (4, 2): Q(1, 2), (6, 0): Q(-1, 4)}
    w2_terms = [("W", W, pa({(1, 0): Q(1)}, ps(Q(-1, 2), L1))),
                ("F", F, {(2, 1): Q(-1)}), ("G5", G5, L5), ("H", H, ps(-1, L6))]

    mixed_residual = {}
    for name, field, time_profile in w2_terms:
        raw = {}
        for p, x in field.items():
            for q, temporal_vec in z3.items():
                k = va(p, q)
                target = raw.setdefault(k, {})
                coefficient = dot(x, q)
                for key, y in temporal_vec.items():
                    # Negative of C(w2,z3): i times +(x dot q)y+(y dot p)x.
                    vpa_into(target, key, va(vs(coefficient, y), vs(dot(y, p), x)))
        for k, temporal_vec in raw.items():
            lam = n2(k)
            if not lam:
                assert all(v == ZERO for v in temporal_vec.values())
                continue
            target = mixed_residual.setdefault(k, {})
            for (r, n), vec in temporal_vec.items():
                projected = va(vec, vs(-dot(k, vec) / lam, k))
                if projected == ZERO:
                    continue
                for (s, m), coefficient in time_profile.items():
                    vpa_into(target, (r + s, n + m), vs(coefficient, projected))
        print(json.dumps({"stage": "mixed_field_constructed", "field": name}), flush=True)
    mixed_residual = {k: {key: v for key, v in profile.items() if v != ZERO}
                      for k, profile in mixed_residual.items()}
    mixed_residual = {k: p for k, p in mixed_residual.items() if p}
    mixed_defect = {k: vector_heat(profile, int(n2(k))) for k, profile in mixed_residual.items()}
    assertions = 0
    for k, profile in mixed_defect.items():
        assert vector_operator(profile, int(n2(k))) == mixed_residual[k]
        assert all(dot(k, value) == 0 for value in profile.values())
        initial = ZERO
        for (r, n), value in profile.items():
            if n == 0:
                initial = va(initial, value)
        assert initial == ZERO
        assertions += 3
    print(json.dumps({"stage": "exact_mixed_heat_complete", "modes": len(mixed_defect),
                      "vector_time_terms": sum(len(p) for p in mixed_defect.values())}), flush=True)

    groups = {}
    mode_records = []
    mixed_D = Q(0)
    for k, profile in sorted(mixed_defect.items()):
        squared = vector_square(profile)
        value = pi(squared)
        assert value >= 0
        upper = sqrt_bounds(value)[1]
        mixed_D += upper
        key = tuple(sorted(squared.items()))
        groups[key] = groups.get(key, 0) + 1
        mode_records.append({"k": list(k), "lambda": qs(n2(k)), "exact_L2_squared": qs(value),
                             "L2_norm_upper": qs(upper),
                             "imaginary_vector_profile": [
                                 {"rate": rate, "power": power, "vector": [qs(x) for x in vec]}
                                 for (rate, power), vec in sorted(profile.items())]})
        assertions += 1
    D3 = Q(previous["global_X_norm_bound"]["rational_upper_bound"])
    c_upper = sqrt_bounds(Q(1, 8))[1]
    quadratic_global = c_upper * D3**2
    full_global = mixed_D + quadratic_global
    print(json.dumps({"stage": "global_bound_complete", "mixed_D_upper": qs(mixed_D),
                      "full_D_upper": qs(round_upper(full_global)),
                      "full_D_decimal": float(full_global), "square_profile_groups": len(groups)}), flush=True)

    grid = [Q(0), Q(1, 16), Q(1, 8), Q(3, 16), Q(1, 4), Q(5, 16), Q(3, 8),
            Q(7, 16), Q(1, 2), Q(5, 8), Q(3, 4), Q(7, 8), Q(1), Q(5, 4),
            Q(3, 2), Q(2), Q(3), None]
    mixed_blocks = [Q(0)] * (len(grid) - 1)
    for key, count in groups.items():
        antiderivative = tail(dict(key))
        endpoints = [endpoint(antiderivative, t) for t in grid]
        for i in range(len(mixed_blocks)):
            integral_upper = endpoints[i][1] - endpoints[i + 1][0]
            assert integral_upper >= 0
            mixed_blocks[i] += count * sqrt_bounds(integral_upper)[1]

    # Pointwise Z>=A0(z3) comes from positive causal field-shell profiles.
    Z = {}
    for term in previous["pointwise_A0_upper_profile"]["terms"]:
        Z = pa(Z, ps(Q(term["amplitude_rational_upper"]), parse_profile(term["causal_heat_profile"])))
    Z2 = pm(Z, Z)
    exact_Z2_integral = pi(Z2)
    assert 0 <= exact_Z2_integral < D3**2
    assertions += 1
    # Integrating the complete pointwise envelope preserves cross-profile
    # cancellations and improves the earlier modewise Minkowski norm bound.
    improved_D3_upper = sqrt_bounds(exact_Z2_integral)[1]
    quadratic_global = c_upper * min(D3**2, exact_Z2_integral)
    full_global = mixed_D + quadratic_global
    z3_spacetime_L2_squared = sum((Q(mode["exact_L2_time_squared"])
                                  for mode in previous["all_defect_modes"]), Q(0))
    z3_actual_lower = sqrt_bounds(z3_spacetime_L2_squared)[0]
    assert round_upper(full_global) < z3_actual_lower
    assertions += 1
    Z2_tail = tail(Z2)
    shifted_tail = tail({(r - 1, n): value for (r, n), value in Z2.items()})
    shifted_total = pi({(r - 1, n): value for (r, n), value in Z2.items()})
    endpoints = [endpoint(Z2_tail, t) for t in grid]
    blocks = []
    for i, (a, b) in enumerate(zip(grid, grid[1:])):
        shifted_at_a = endpoint(shifted_tail, a)
        past_integral_upper = max(Q(0), shifted_total - shifted_at_a[0])
        weighted_past_upper = expneg(a)[1] * past_integral_upper
        inside_upper = max(Q(0), endpoints[i][1] - endpoints[i + 1][0])
        history_quad = c_upper * (weighted_past_upper + inside_upper)
        quadratic_upper = min(history_quad, quadratic_global)
        blocks.append({"left": qs(a), "right": "infinity" if b is None else qs(b),
                       "mixed_defect_norm_upper": qs(mixed_blocks[i]),
                       "quadratic_defect_norm_upper": qs(round_upper(quadratic_upper)),
                       "full_defect_norm_upper": qs(round_upper(mixed_blocks[i] + quadratic_upper)),
                       "full_upper_decimal_for_orientation": float(mixed_blocks[i] + quadratic_upper)})

    result = {
        "status": "PASS_FOURTH_HEAT_DEFECT_STATIC_CHECKS",
        "scope": "Exact mixed Fourier/exponential-polynomial construction plus rigorous bilinear quadratic envelopes. No ODE/PDE integration or numerical quadrature. Historical z3 JSON read as data only; no old scripts imported or executed.",
        "source": {"filename": source.name, "sha256": hashlib.sha256(source_bytes).hexdigest()},
        "definitions": {"r3": "-C(w2,z3)+B(z3,z3)", "z4": "H r3", "mixed": "-H C(w2,z3)",
                        "quadratic": "H B(z3,z3)", "w3": "w2-z3", "amplitude_a": "1/1"},
        "mixed_mode_count": len(mixed_defect),
        "max_mixed_squared_wavenumber": qs(max(n2(k) for k in mixed_defect)),
        "global_bounds": {
            "mixed_X_upper": qs(mixed_D), "mixed_upper_decimal": float(mixed_D),
            "D3_source_upper": qs(D3), "bilinear_c_upper": qs(c_upper),
            "pointwise_Z_exact_integral_squared": qs(exact_Z2_integral),
            "pointwise_Z_integral_squared_decimal": float(exact_Z2_integral),
            "improved_D3_X_upper": qs(improved_D3_upper),
            "quadratic_X_upper": qs(round_upper(quadratic_global)),
            "full_z4_X_upper": qs(round_upper(full_global)), "full_upper_decimal": float(full_global),
            "mixed_bound_method": "Minkowski sum of exact modewise L2 norms, with rational square-root upper bounds.",
            "quadratic_bound_method": "||H B(z3,z3)||X <= (1/(2sqrt2)) ||z3||X^2 <= c_upper min(D3_source_upper^2, integral Z^2).",
        },
        "actual_single_step_decrease": {
            "z3_spacetime_L2_exact_squared": qs(z3_spacetime_L2_squared),
            "z3_X_rational_lower_bound": qs(z3_actual_lower),
            "z4_X_rational_upper_bound": qs(round_upper(full_global)),
            "ratio_rational_upper_bound": qs(round_upper(full_global) / z3_actual_lower),
            "argument": "A0(z3(t)) >= (sum_k |z3_k(t)|^2)^(1/2), so ||z3||X >= space-time L2 norm. The resulting lower bound exceeds the rigorous z4 upper bound.",
            "scope": "A single actual decrease ||z4||X < ||z3||X is proved; this does not establish a contraction for later iterates or convergence of the whole iteration.",
        },
        "block_defect_bounds": blocks,
        "quadratic_block_method": "On [a,b], norm <= c[int_0^a exp(-(a-s))Z(s)^2 ds+int_a^b Z(s)^2 ds], where pointwise Z>=A0(z3). Take minimum with global bound. Full past history retained.",
        "mixed_mode_profiles": mode_records,
        "pointwise_Z_upper": pjson(Z),
        "number_of_assertions": assertions,
        "qualification": "Only one additional Picard defect is bounded. No conclusion about convergence of all later iterates follows from a finite decrease or these bounds alone.",
    }
    destination = Path(__file__).with_name("fourth_heat_defect_static_result.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "assertions": assertions,
                      "global_bounds": result["global_bounds"],
                      "block_upper_decimals": [b["full_upper_decimal_for_orientation"] for b in blocks],
                      "result_file": str(destination)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
