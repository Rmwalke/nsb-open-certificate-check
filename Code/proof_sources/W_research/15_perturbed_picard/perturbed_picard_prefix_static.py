"""Fresh rational finite-prefix entry check; no evolution or quadrature.

Read frozen certificates as data. Derive the heat-seeded w0h,w1h,w2h,w3h
differences from exact endpoint moments and causal bilinear block bounds.
No historical checker is imported or executed.
"""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json

SOURCE = Path("perturbation_block_gate_result.json")
BASE = Path("refined_W_block_gate_result.json")
source_bytes = SOURCE.read_bytes()
base_bytes = BASE.read_bytes()
data = json.loads(source_bytes)
base = json.loads(base_bytes)
assert data["source_sha256"] == hashlib.sha256(base_bytes).hexdigest()
certificate = data["selected_certificate"]
assert certificate["passed"]
unit = [t for t in base["tests"]
        if t["name"] == "w3_prior17_exact_U2_plus_D3"][0]
assert unit["passed"]

SCALE = 10**30
ESCALE = 10**90
checks = 0

def check(condition):
    global checks
    assert condition
    checks += 1

def lower(q, scale):
    return F(q.numerator*scale//q.denominator, scale)

def upper(q, scale=SCALE):
    return F((q.numerator*scale+q.denominator-1)//q.denominator, scale)

def exp_minus(q):
    check(q >= 0)
    squarings = 0
    while q > 1:
        q /= 2
        squarings += 1
    even = sum(((-q)**n/F(factorial(n)) for n in range(81)), F(0))
    odd = even- q**81/F(factorial(81))
    lo, hi = lower(odd, ESCALE), upper(even, ESCALE)
    check(0 < lo <= hi)
    for _ in range(squarings):
        lo, hi = lower(lo*lo, ESCALE), upper(hi*hi, ESCALE)
    return lo, hi

def sqrt_upper(q):
    check(q >= 0)
    floor_root = isqrt(q.numerator*SCALE*SCALE//q.denominator)
    cap = F(floor_root, SCALE)
    if cap*cap < q:
        cap += F(1, SCALE)
    check(cap*cap >= q)
    return cap

def reference_tail(t, order):
    if t is None:
        return F(0), F(0)
    terms = [(F(9,2), 2)]
    if order == 1:
        terms += [(6*t+2, 3),
                  (F(9,4)*t*t+F(9,8)*t+F(9,32), 4)]
    lo = hi = F(0)
    for coefficient, rate in terms:
        el, eh = exp_minus(rate*t)
        lo += coefficient*el
        hi += coefficient*eh
    return lo, hi

grid = [None if t == "infinity" else F(t) for t in data["grid"]]
check(grid == [None if t == "infinity" else F(t) for t in unit["grid"]])
n = len(grid)-1
check(n == 17 and grid[0] == 0 and grid[-1] is None)
c = F(data["c_upper"])
chi = F(data["selected_chi"])
check(chi == F(1,10000))
check(c*c > F(1,8))
K = [[F(x) for x in row] for row in data["source_memory_factors"]]
for i in range(n):
    check(len(K[i]) == i+1 and K[i][-1] == 1)
    check(K[i][:-1] == [F(x) for x in unit["rows"][i]["coefficients"]])

# Exact moment caps for the first three unperturbed Picard iterates.
# w0=e^-t W; A0(w1)=3e^-t+3t e^-2t; V2 comes from frozen exact U2.
V = [[], [], []]
for i, (a,b) in enumerate(zip(grid[:-1], grid[1:])):
    for order in (0,1):
        left = reference_tail(a, order)
        right = reference_tail(b, order)
        V[order].append(sqrt_upper(left[1]-right[0]))
    V[2].append(F(unit["rows"][i]["reference2_norm_upper"]))

s = [chi*F(x) for x in data["heat_profile_caps"]]
rho = [F(x) for x in certificate["radii_upper"]]
check(s[0]**2 >= chi*chi/2)
for i in range(n):
    check(s[i] == F(certificate["rows"][i]["initial_heat_cap"]))

# E[k][i] bounds ||w_k,h-w_k|| on block i, with all source history retained.
E = [s]
for order in range(3):
    current = []
    for i in range(n):
        bound = s[i]+c*sum((K[i][j]*(2*V[order][j]*E[order][j]
                              +E[order][j]**2)
                            for j in range(i+1)), F(0))
        cap = upper(bound)
        check(cap >= bound)
        current.append(cap)
    E.append(current)

margins = [rho[i]-E[3][i] for i in range(n)]
for margin in margins:
    check(margin > 0)
entry_ratios = [E[3][i]/rho[i] for i in range(n)]
max_ratio = max(entry_ratios)
check(max_ratio < 1)

matrix = [[F(x) for x in row] for row in certificate["contraction_matrix"]]
q = F(3,4)
omega = []
for i in range(n):
    check(matrix[i][i] < q)
    omega.append((1+sum((matrix[i][j]*omega[j] for j in range(i)), F(0)))
                 /(q-matrix[i][i]))
for i in range(n):
    check(sum((matrix[i][j]*omega[j] for j in range(i+1)), F(0))
          == q*omega[i]-1)

# Initial fixed-point error is bounded by E3+rho, not merely rho.
initial_error = [E[3][i]+rho[i] for i in range(n)]
weighted_initial = max(initial_error[i]/omega[i] for i in range(n))
rate_constant_squared = weighted_initial**2*sum((x*x for x in omega), F(0))
rate_constant = sqrt_upper(rate_constant_squared)

result = {
    "status": "PASS_PERTURBED_TRUE_HEAT_SEED_FINITE_PREFIX_ENTRY",
    "scope": "Smooth real mean-zero divergence-free h; A_-1(h)<=1/10000; full initial datum W+h.",
    "source_json": str(SOURCE),
    "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
    "baseline_source_json": str(BASE),
    "baseline_source_sha256": hashlib.sha256(base_bytes).hexdigest(),
    "assertions_passed": checks,
    "chi_radius": str(chi),
    "grid": data["grid"],
    "unperturbed_reference_caps": [list(map(str,row)) for row in V],
    "prefix_difference_caps": [list(map(str,row)) for row in E],
    "entry_margins": list(map(str,margins)),
    "maximum_entry_ratio": str(max_ratio),
    "maximum_entry_ratio_for_orientation": float(max_ratio),
    "initial_fixed_point_error_caps": list(map(str,initial_error)),
    "rate": str(q),
    "global_rate_prefactor_upper": str(rate_constant),
    "global_rate_prefactor_for_orientation": float(rate_constant),
    "error_statement": "For the true heat-seeded Picard sequence v[n], n>=3: block error <=M^(n-3)(E3+rho), and global X error <=C*(3/4)^(n-3).",
    "not_claimed": ["Enlarged global regularity class", "A_-1(h)<=1/20 Picard convergence", "All smooth data convergence"]
}
Path("perturbed_picard_prefix_static_result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:result[k] for k in ["status","assertions_passed","chi_radius",
    "maximum_entry_ratio_for_orientation","global_rate_prefactor_for_orientation"]},indent=2))
