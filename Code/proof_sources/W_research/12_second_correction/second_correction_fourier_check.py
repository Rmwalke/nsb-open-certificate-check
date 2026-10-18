"""Fresh exact Fourier algebra for a second nonlinear heat correction.

No old files/code are loaded, no evolution is simulated, and no quadrature is
performed. Rational coefficients represent the imaginary part of each mode.
"""

from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json


ZERO = (Q(0), Q(0), Q(0))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def norm2(a):
    return dot(a, a)


def qs(q):
    q = Q(q)
    return f"{q.numerator}/{q.denominator}"


def field_scale(c, field):
    return {k: scale(c, value) for k, value in field.items() if c}


def field_add(u, v):
    result = dict(u)
    for k, value in v.items():
        result[k] = add(result.get(k, ZERO), value)
    return {k: value for k, value in result.items() if value != ZERO}


def B(u, v):
    raw = {}
    for p, up in u.items():
        for q, vq in v.items():
            k = add(p, q)
            raw[k] = add(raw.get(k, ZERO), scale(-dot(up, q), vq))
    result = {}
    for k, value in raw.items():
        k2 = norm2(k)
        if not k2:
            assert value == ZERO
            continue
        projected = add(value, scale(-dot(k, value) / k2, k))
        if projected != ZERO:
            result[k] = projected
    return result


def C(u, v):
    return field_add(B(u, v), B(v, u))


def radical(q):
    q = Q(q)
    if not q:
        return Q(0), 1
    n, outside, inside, p = q.numerator * q.denominator, 1, 1, 2
    while p * p <= n:
        power = 0
        while n % p == 0:
            n //= p
            power += 1
        outside *= p ** (power // 2)
        inside *= p ** (power % 2)
        p += 1
    if n > 1:
        inside *= n
    coefficient = Q(outside, q.denominator)
    assert coefficient**2 * inside == q
    return coefficient, inside


def norm_record(field, critical):
    terms = {}
    for k, value in field.items():
        c, d = radical(norm2(value) / norm2(k) if critical else norm2(value))
        terms[d] = terms.get(d, Q(0)) + c
    denominator = 10**9
    lower = upper = Q(0)
    for d, c in terms.items():
        n = isqrt(d * denominator**2)
        lo = Q(n, denominator)
        hi = lo if n**2 == d * denominator**2 else Q(n + 1, denominator)
        assert lo**2 <= d <= hi**2
        lower += c * lo
        upper += c * hi
    return {
        "exact_radical_sum": [{"coefficient": qs(c), "sqrt_integer": d} for d, c in sorted(terms.items())],
        "rational_lower_bound": qs(lower),
        "rational_upper_bound": qs(upper),
        "upper_decimal_for_orientation_only": float(upper),
    }


def field_record(field):
    counts = {}
    for k, value in field.items():
        key = (norm2(k), norm2(value))
        counts[key] = counts.get(key, 0) + 1
    return {
        "mode_count": len(field),
        "shell_counts": [{"k_squared": qs(k2), "coefficient_norm_squared": qs(v2), "count": n}
                         for (k2, v2), n in sorted(counts.items())],
        "critical_norm": norm_record(field, True),
        "A0_norm": norm_record(field, False),
        "fundamental_modes": [{"k": list(k), "imaginary_rational_vector": [qs(x) for x in value]}
                              for k, value in sorted(field.items()) if norm2(k) == 1],
        "all_coefficients": [{"k": list(k), "imaginary_rational_vector": [qs(x) for x in value]}
                             for k, value in sorted(field.items())],
    }


def linear_heat_operator(profile, eigenvalue):
    """Apply d/dt+eigenvalue to a finite sum polynomial(t)*exp(-rate*t)."""
    result = {}
    for rate, coefficients in profile.items():
        out = [(eigenvalue - rate) * coefficient for coefficient in coefficients]
        for j in range(1, len(coefficients)):
            out[j - 1] += j * coefficients[j]
        while out and out[-1] == 0:
            out.pop()
        if out:
            result[rate] = tuple(out)
    return result


def main():
    checks = []

    def check(name, result):
        assert result, name
        checks.append(name)

    W = {}
    for wave_axis, component in ((2, 0), (0, 1), (1, 2)):
        for sign in (-1, 1):
            k = tuple(sign if i == wave_axis else 0 for i in range(3))
            W[k] = tuple(Q(-sign, 2) if i == component else Q(0) for i in range(3))
    F = B(W, W)
    G = C(W, F)
    G1 = {k: value for k, value in G.items() if norm2(k) == 1}
    G5 = {k: value for k, value in G.items() if norm2(k) == 5}
    H = B(F, F)
    check("G_is_G1_plus_G5", G == field_add(G1, G5))
    check("G1_is_minus_W_over_two", G1 == field_scale(Q(-1, 2), W))
    check("H_shell_six", len(H) == 24 and all(norm2(k) == 6 for k in H))
    check("W_F_G5_H_have_pairwise_disjoint_shells", all(
        set(x).isdisjoint(y)
        for i, x in enumerate((W, F, G5, H)) for y in (W, F, G5, H)[i + 1:]
    ))

    independent = {"W": W, "F": F, "G5": G5, "H": H}
    pairs = {}
    keys = list(independent)
    for i, left in enumerate(keys):
        for right in keys[i:]:
            key = left + "," + right
            pairs[key] = B(independent[left], independent[right]) if left == right else C(independent[left], independent[right])
    for name, field in {**independent, "G1": G1, "G": G, **pairs}.items():
        check(name + "_mean_zero", (0, 0, 0) not in field)
        check(name + "_solenoidal", all(dot(k, value) == 0 for k, value in field.items()))
        check(name + "_real", all(field[tuple(-x for x in k)] == scale(-1, value) for k, value in field.items()))

    coarse_caps = {"W,W": Q(17, 8), "W,F": Q(87, 40), "W,G5": Q(171, 100),
                   "W,H": Q(27, 20), "F,F": Q(71, 100), "F,G5": Q(6, 5),
                   "F,H": Q(183, 100), "G5,G5": Q(11, 25), "G5,H": Q(3, 4),
                   "H,H": Q(9, 25)}
    for name, cap in coarse_caps.items():
        check(name + "_coarse_rational_critical_cap",
              Q(norm_record(pairs[name], True)["rational_upper_bound"]) < cap)
    low_G5H = {k: value for k, value in pairs["G5,H"].items() if norm2(k) == 1}
    check("C_G5_H_fundamental_part_W_over_24", low_G5H == field_scale(Q(1, 24), W))
    check("no_other_pair_fundamental_modes", all(
        all(norm2(k) != 1 for k in field) for name, field in pairs.items()
        if name not in ("W,F", "G5,H")
    ))

    L1 = {1: (Q(1, 4),), 3: (Q(-1, 4), Q(-1, 2))}
    L5 = {5: (Q(1, 4),), 3: (Q(-1, 4), Q(1, 2))}
    L6 = {4: (Q(1, 4), Q(-1, 2), Q(1, 2)), 6: (Q(-1, 4),)}
    for name, profile, eigenvalue, expected in (
        ("L1", L1, 1, {3: (Q(0), Q(1))}),
        ("L5", L5, 5, {3: (Q(0), Q(1))}),
        ("L6", L6, 6, {4: (Q(0), Q(0), Q(1))}),
    ):
        check(name + "_zero_initial_value", sum(poly[0] for poly in profile.values()) == 0)
        check(name + "_exact_forced_heat_identity", linear_heat_operator(profile, eigenvalue) == expected)

    # The following algebra is checked at independent rational formal samples.
    # It supplements the displayed exact bilinear expansion; it is not an
    # assertion that finitely many samples prove a polynomial identity.
    for alpha, beta, gamma, eta, theta in ((Q(1), Q(-2), Q(-3), Q(4), Q(-5)),
                                         (Q(2, 3), Q(-5, 7), Q(-11, 13), Q(17, 19), Q(-23, 29))):
        w1 = field_add(field_scale(alpha, W), field_scale(beta, F))
        correction = field_add(field_add(field_scale(gamma, W), field_scale(eta, G5)), field_scale(theta, H))
        direct = field_add(C(w1, correction), B(correction, correction))
        terms = (
            (2 * alpha * gamma + gamma**2, pairs["W,W"]),
            (beta * gamma, pairs["W,F"]),
            ((alpha + gamma) * eta, pairs["W,G5"]),
            (beta * eta, pairs["F,G5"]),
            ((alpha + gamma) * theta, pairs["W,H"]),
            (beta * theta, pairs["F,H"]),
            (eta**2, pairs["G5,G5"]),
            (theta**2, pairs["H,H"]),
            (eta * theta, pairs["G5,H"]),
        )
        expanded = {}
        for coefficient, field in terms:
            expanded = field_add(expanded, field_scale(coefficient, field))
        check("direct_residual_matches_expansion_at_" + qs(alpha), direct == expanded)

    result = {
        "status": "PASS_SECOND_CORRECTION_FOURIER_STATIC_CHECKS",
        "scope": "Fresh self-contained exact Fourier arithmetic; no old data/code, time integration, quadrature, formal proof or unknown-flow numerical certificate.",
        "conventions": {
            "B": "P[(u dot grad)v]",
            "distinct_pair": "C(A,B)=B(A,B)+B(B,A)",
            "diagonal_pair": "B(A,A), with no factor two",
            "G1": "-W/2",
            "G5": "G-G1",
        },
        "profiles": {
            "L1": "[exp(-t)-(1+2t)exp(-3t)]/4",
            "L5": "[exp(-5t)+(2t-1)exp(-3t)]/4",
            "L6": "[(2t^2-2t+1)exp(-4t)-exp(-6t)]/4",
            "all_L_nonnegative": "Each is a causal heat convolution of a nonnegative scalar forcing.",
            "alpha": "a exp(-t)",
            "beta": "-a^2 t exp(-2t)",
            "gamma": "-a^3 L1/2",
            "eta": "a^3 L5",
            "theta": "-a^4 L6",
            "w1": "alpha W+beta F",
            "correction": "gamma W+eta G5+theta H",
            "w2": "(alpha+gamma)W+beta F+eta G5+theta H",
            "z2_parent_sign_convention": "z2=w1-w2=-correction, so r2=-C(w1,z2)+B(z2,z2).",
        },
        "base_fields": {key: field_record(value) for key, value in {**independent, "G1": G1, "G": G}.items()},
        "pair_fields": {key: field_record(value) for key, value in pairs.items()},
        "coarse_rational_pair_critical_upper_bounds": {key: qs(value) for key, value in coarse_caps.items()},
        "residual_expansion": [
            {"coefficient": "2alpha gamma+gamma^2", "pair": "W,W"},
            {"coefficient": "beta gamma", "pair": "W,F"},
            {"coefficient": "(alpha+gamma)eta", "pair": "W,G5"},
            {"coefficient": "beta eta", "pair": "F,G5"},
            {"coefficient": "(alpha+gamma)theta", "pair": "W,H"},
            {"coefficient": "beta theta", "pair": "F,H"},
            {"coefficient": "eta^2", "pair": "G5,G5"},
            {"coefficient": "theta^2", "pair": "H,H"},
            {"coefficient": "eta theta", "pair": "G5,H"},
        ],
        "triangle_majorant": "R2(t)<=sum over displayed residual_expansion of |coefficient(t)| chi(pair_field)",
        "positive_profile_majorants": {
            "definitions": "p=a exp(-t), q=a^2 t exp(-2t), d=a^3 L1/2, e=a^3 L5, f=a^4 L6; all nonnegative for a>=0",
            "range": "0<=a<=sqrt(3), so d<=p via L1<=exp(-t)/4",
            "retaining_cancellations": "chi(F)(2pd-d^2)+chi(G)qd+chi(C(W,G5))(p-d)e+chi(C(F,G5))qe+chi(C(W,H))(p-d)f+chi(C(F,H))qf+chi(B(G5,G5))e^2+chi(B(H,H))f^2+chi(C(G5,H))ef",
            "fully_positive_polynomial_upper_bound": "2chi(F)pd+chi(G)qd+chi(C(W,G5))pe+chi(C(F,G5))qe+chi(C(W,H))pf+chi(C(F,H))qf+chi(B(G5,G5))e^2+chi(B(H,H))f^2+chi(C(G5,H))ef",
            "comparison": "The fully positive bound drops only the nonpositive terms -chi(F)d^2, -chi(C(W,G5))de and -chi(C(W,H))df.",
        },
        "exact_fundamental_residual": "r2_fundamental=-(qd/2+ef/24)W; G_fundamental=-W/2 and C(G5,H)_fundamental=W/24, with both contributions having the same sign.",
        "exact_reference_amplitude_when_0_le_a_le_sqrt8": "A0(w2)=3(alpha+gamma)-3beta+(3/2)eta-sqrt(3)theta; alpha+gamma>=a exp(-t)(1-a^2/8)>=0",
        "number_of_assertions": len(checks),
        "assertions": checks,
    }
    destination = Path(__file__).with_name("second_correction_fourier_check_result.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "number_of_assertions": len(checks),
        "pairs": {key: {"mode_count": len(value), "chi": norm_record(value, True),
                        "fundamental_mode_count": sum(norm2(k) == 1 for k in value)} for key, value in pairs.items()},
        "result_file": str(destination),
    }, indent=2))


if __name__ == "__main__":
    main()
