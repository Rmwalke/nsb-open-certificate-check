"""Fresh exact scalar algebra for the second corrected W reference.

No historical research module is imported or executed. No ODE/PDE is
integrated. All time integrals below are exact exponential-polynomial
factorial identities, followed by finite positive-series tail bounds.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

checks = []

def check(name, condition):
    assert condition, name
    checks.append(name)

def add(*polys):
    result = {}
    for poly in polys:
        for key, value in poly.items():
            result[key] = result.get(key, Q(0)) + value
    return {key: value for key, value in result.items() if value}

def scale(poly, value):
    return {key: coefficient * value for key, coefficient in poly.items()
            if coefficient * value}

def multiply(first, second):
    result = {}
    for (rate, degree), coefficient in first.items():
        for (other_rate, other_degree), other in second.items():
            key = (rate + other_rate, degree + other_degree)
            result[key] = result.get(key, Q(0)) + coefficient * other
    return {key: value for key, value in result.items() if value}

def moment(poly, series_index=0):
    return sum((coefficient * factorial(degree) /
                Q((rate + 2 * series_index) ** (degree + 1))
                for (rate, degree), coefficient in poly.items()), Q(0))

def exponential_upper(value):
    term = result = Q(1)
    for index in range(1, 17):
        term = term * value / index
        result += term
    check('Taylor tail ratio below one: ' + str(value), value < 18)
    return result + (term * value / 17) / (1 - value / 18)

L1 = {(1, 0): Q(1, 4), (3, 0): Q(-1, 4), (3, 1): Q(-1, 2)}
L5 = {(5, 0): Q(1, 4), (3, 0): Q(-1, 4), (3, 1): Q(1, 2)}
L6 = {(4, 2): Q(1, 2), (4, 1): Q(-1, 2),
      (4, 0): Q(1, 4), (6, 0): Q(-1, 4)}

field_names = ['F', 'K_WC', 'K_WH', 'G', 'K_FC', 'K_FH',
               'B_CC', 'K_CH', 'B_HH']
caps = [Q(17, 8), Q(171, 100), Q(27, 20), Q(87, 40), Q(6, 5),
        Q(183, 100), Q(11, 25), Q(3, 4), Q(9, 25)]
# Certified rational bounds supplied by the separate new exact Fourier check.
fourier_bounds = [Q(4242640689, 2000000000),
                  Q(34058305149, 20000000000),
                  Q(442525417063, 330000000000),
                  Q(10854101967, 5000000000),
                  Q(52273526363, 44000000000),
                  Q(764598478601, 420000000000),
                  Q(12059252349, 28000000000),
                  Q(9660978826560959, 12932920000000000),
                  Q(543918125957, 1540000000000)]
for name, certified, cap in zip(field_names, fourier_bounds, caps):
    check('Fourier cap ' + name, certified < cap)

def certificate(amplitude):
    a = amplitude
    x = {(1, 0): a}
    y = {(2, 1): a * a}
    p = scale(L1, a ** 3 / 2)
    q = scale(L5, a ** 3)
    r = scale(L6, a ** 4)
    x_minus_p = add(x, scale(p, -1))
    profiles = [add(scale(multiply(x, p), 2), scale(multiply(p, p), -1)),
                multiply(x_minus_p, q), multiply(x_minus_p, r),
                multiply(y, p), multiply(y, q), multiply(y, r),
                multiply(q, q), multiply(q, r), multiply(r, r)]
    majorant = add(*(scale(profile, cap)
                    for profile, cap in zip(profiles, caps)))
    c = Q(9, 8) * a * a
    delta = a ** 3 / 2 + Q(9, 128) * a ** 4
    reference_budget_upper = Q(9, 2) * a * a + 2 * a ** 3 + Q(9, 32) * a ** 4
    multiplier = Q(1)
    terms = []
    for index in range(10):
        if index:
            multiplier = multiplier * c / index
        current_moment = moment(majorant, index)
        check('positive exact moment %s/%s' % (a, index), current_moment > 0)
        terms.append(multiplier * current_moment)
    check('series remainder ratio below one ' + str(a), c < 10)
    series_upper = sum(terms[:9], Q(0)) + terms[9] / (1 - c / 10)
    return {'amplitude': a, 'total_residual_upper': moment(majorant),
            'series_upper': series_upper, 'delta': delta,
            'exp_delta_upper': exponential_upper(delta),
            'reference_budget_upper': reference_budget_upper,
            'initial_multiplier_upper': exponential_upper(reference_budget_upper / 4),
            'residual_endpoint_upper': series_upper * exponential_upper(delta)}

A = certificate(Q(11, 10))
B = certificate(Q(1))
check('segment lies below sqrt3', Q(11, 10) ** 2 < 3)
check('monotonicity derivative F', 4 - Q(6, 16) * 3 > 0)
check('monotonicity derivative WC', 4 - Q(6, 8) * 3 > 0)
check('monotonicity derivative WH', 5 - Q(7, 8) * 3 > 0)
check('A residual integral <3/20', A['total_residual_upper'] < Q(3, 20))
check('A series <27/125', A['series_upper'] < Q(27, 125))
check('A expDelta <54/25', A['exp_delta_upper'] < Q(54, 25))
check('A endpoint residual <47/100', Q(27, 125) * Q(54, 25) < Q(47, 100))
check('A initial multiplier <17/2', A['initial_multiplier_upper'] < Q(17, 2))
check('A P <213/25', A['reference_budget_upper'] < Q(213, 25))
check('A initial radius gate', Q(47, 100) + Q(17, 2) / 100 == Q(111, 200))
check('A gate <3/5', Q(111, 200) < Q(3, 5))
check('B residual integral <101/1000', B['total_residual_upper'] < Q(101, 1000))
check('B series <271/2000', B['series_upper'] < Q(271, 2000))
check('B expDelta <177/100', B['exp_delta_upper'] < Q(177, 100))
check('B endpoint residual <6/25', Q(271, 2000) * Q(177, 100) < Q(6, 25))
check('B initial multiplier <11/2', B['initial_multiplier_upper'] < Q(11, 2))
check('B initial radius gate', Q(6, 25) + Q(11, 2) / 20 == Q(103, 200))
check('B gate <3/5', Q(103, 200) < Q(3, 5))
check('Psi(1)>3/5', exponential_upper(Q(1, 2)) < Q(5, 3))
check('A F <4/25', Q(1, 100) + Q(3, 20) == Q(4, 25))
check('B F <4/25', Q(1, 20) + Q(101, 1000) < Q(4, 25))
check('H<35', 4 * Q(213, 25) + 4 * Q(4, 25) == Q(868, 25) < 35)
check('P<9', Q(213, 25) < 9)
check('sqrtH<6', 35 < 36)
check('S2 through sqrt6', 6 + Q(12, 2) + Q(8, 3) + Q(6, 4) +
      Q(24, 5) + Q(24, 6) == Q(749, 30) < 25)
check('sqrt(3/2)<5/4', Q(3, 2) < Q(5, 4) ** 2)
check('sqrt6<5/2', 6 < Q(5, 2) ** 2)
check('A low projection norm <7/5', Q(11, 10) * Q(5, 4) + Q(5, 2) / 100 == Q(7, 5))
check('B low projection norm <7/5', Q(5, 4) + Q(5, 2) / 20 < Q(7, 5))
check('low budget <81', Q(25, 2) * Q(7, 5) ** 2 == Q(49, 2) < 81)

def encode(value):
    if isinstance(value, Q):
        return {'numerator': value.numerator, 'denominator': value.denominator,
                'decimal_display_only': float(value)}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value

result = {'status': 'PASS_SECOND_CORRECTED_W_ANALYTIC_GATE',
          'assertion_count': len(checks), 'checks': checks,
          'scope': 'Exact factorial moments and positive-series tails; no ODE/PDE evolution.',
          'case_A': A, 'case_B': B,
          'proved_classes': ['|a|<=11/10, A_-1(h)<=1/100',
                            '|a|<=1, A_-1(h)<=1/20'],
          'uniform_error_barrier': Q(1),
          'integrated_A1_error_bound': Q(35),
          'all_cutoff_A0_squared_integral_bound': Q(81),
          'NSB_bound': '(E0*X0/2)*exp(81)'}
Path(__file__).with_name('second_corrected_gate_static_result.json').write_text(
    json.dumps(encode(result), indent=2) + '\n')
print(result['status'], len(checks), 'exact assertions')
for label, data in [('A', A), ('B', B)]:
    print(label, 'residual_endpoint_upper', float(data['residual_endpoint_upper']),
          'initial_multiplier_upper', float(data['initial_multiplier_upper']))
