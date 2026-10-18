"""Recheck displayed rational margins; no ODE integration and no network."""
import hashlib
import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

BASE = Path(__file__).resolve().parent
r = json.loads((BASE / 'n3_validated_result.json').read_text())
assert r['status'] == 'PASS_ORIGINAL_GATE'
assert r['completed_horizon'] == '8'
assert hashlib.sha256((BASE / 'n3_exact_algebra_intermediate.json').read_bytes()).hexdigest() == r['source_sha256']
last = r['records'][-1]
Q = 1 << r['precision_bits']
eta = F(last['eta_prefix'][1], Q) + F(r['eta_tail_upper'])
integ = F(last['I_prefix'][1], Q) + F(r['I_tail_upper'])
factor = 27 * F(25, 3)**6
energy = F(r['terminal_energy_upper'])
den = 10**30
scaled = 660*energy*den**2
root = isqrt(scaled.numerator // scaled.denominator)
if F(root**2) < scaled:
    root += 1
H = F(root, den)
lambda_infty = F(last['lambda'][1], Q) + 2*H + H**2

# Positive rational exponential series with an explicit geometric tail.
x = F(176, 25)
term = total = F(1)
for n in range(1, 129):
    term *= x/n
    total += term
exp_upper = total + (term*x/129)/(1-x/130)
rounded_eta = F(587, 10**9)
rounded_I = F(156400)
radius = F(1, 2500)
checks = {
    'eta_lt_587_over_1e9': eta < rounded_eta,
    'I_lt_156400': integ < rounded_I,
    'recorded_kappa_bounds_exact_product': factor*eta**2*integ <= F(r['kappa_original_upper']),
    'kappa_lt_49_over_100': F(r['kappa_original_upper']) < F(49, 100),
    'rounded_gate_lt_49_over_100': factor*rounded_eta**2*rounded_I < F(49, 100),
    'Lambda_infty_lt_176_over_25': lambda_infty < F(176, 25),
    'exp_176_over_25_lt_1150': exp_upper < 1150,
    'sqrt_1_minus_kappa_gt_7_over_10': 1-F(49,100) > F(7,10)**2,
    'error_constant_lt_1_over_1000': rounded_eta*1150/F(7,10) < F(1,1000),
    'uniform_H1_squared_lt_14': F(27,2) < F(37,10)**2 and F(1,1000) < F(1,30)**2 and (F(37,10)+F(1,30))**2 < 14,
    'neighborhood_radius_1_over_2500_gate_lt_4_over_5': factor*(rounded_eta+radius**2)**2*rounded_I < F(4,5),
}
assert all(checks.values()), checks
out = {'method': 'exact rational inequalities and separate positive-series exponential bound',
       'checks': checks, 'eta_upper': str(eta), 'I_upper': str(integ),
       'Lambda_infty_upper': str(lambda_infty), 'rounded_eta': str(rounded_eta),
       'rounded_I': str(rounded_I), 'kappa_strict_bound': '49/100',
       'H1_neighborhood_radius': str(radius)}
(BASE / 'n3_rational_margin_check.json').write_text(json.dumps(out, indent=2))
print(json.dumps({'all_checks_pass': True, 'checks': checks}, indent=2))
