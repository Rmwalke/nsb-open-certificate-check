"""Fresh rational consequences of passing endpoint gates.

No old scripts are imported or executed. All logarithms and radicals
are enclosed by finite rational expressions; no quadrature/evolution.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
SOURCE = Path("enlarged_endpoint_gate_result.json")
source_bytes = SOURCE.read_bytes()
EXPECTED_SHA = "661a9c1c426b59e5f9d0a9477628f03f3b96a1ec8f86dc69f4c7d5e919ae4f0f"
data = json.loads(source_bytes)
checks = []
SCALE = 10**30

def check(condition, description):
    assert condition, description
    checks.append(description)

check(hashlib.sha256(source_bytes).hexdigest() == EXPECTED_SHA,
      "frozen endpoint source hash")

def atanh_log_interval(x):
    # x in [1,2], log(x)=2 sum q^(2n+1)/(2n+1), 0<=q<=1/3.
    check(F(1) <= x <= 2, "normalized logarithm argument")
    q = (x-1)/(x+1)
    count = 40
    lo = 2*sum((q**(2*n+1)/F(2*n+1) for n in range(count)), F(0))
    hi = lo + 2*q**(2*count+1)/(F(2*count+1)*(1-q*q))
    check(0 <= lo <= hi, "positive logarithm tail enclosure")
    return lo, hi

LOG2 = atanh_log_interval(F(2))

def log_interval(x):
    check(x >= 1, "log argument at least one")
    k = 0
    while x > 2:
        x /= 2
        k += 1
    lo, hi = atanh_log_interval(x)
    return lo+k*LOG2[0], hi+k*LOG2[1]

def primitive_interval(s):
    check(0 <= s < 2, "critical error within increasing endpoint branch")
    lo, hi = log_interval(1/(1-s/2))
    result = (-2*s+4*lo, -2*s+4*hi)
    check(result[1] >= 0, "primitive upper bound nonnegative")
    return result

def sqrt_upper(q):
    check(q >= 0, "nonnegative radical")
    n = isqrt(q.numerator*SCALE*SCALE//q.denominator)
    r = F(n,SCALE)
    if r*r < q:
        r += F(1,SCALE)
    check(r*r >= q, "outward radical enclosure")
    return r

def compute_budget(label, amplitude, radius, barrier, P, R):
    delta = 1-barrier/2
    fb = primitive_interval(barrier)[1]
    fr = primitive_interval(radius)[1]
    C = fr+barrier*R/delta
    rootQ = sqrt_upper(fb)+sqrt_upper(fb+C)
    budget = (sqrt_upper(P)+rootQ)**2
    exponent = (budget.numerator+budget.denominator-1)//budget.denominator
    check(budget < exponent, "strict reported integer budget: "+label)
    check(exponent >= 2, "integer exponent at least two: "+label)
    # e^2 > 1+2+2+8/6=19/3>6; low cutoffs |k|<sqrt6 fit.
    check(F(19,3) > 6, "all lower cutoffs absorbed: "+label)
    # Compare the strongest constant-barrier quadratic estimate as a diagnostic.
    rootH = (sqrt_upper(barrier*P)
             +sqrt_upper(barrier*P+4*delta*(radius+R)))/(2*delta)
    old_budget = (sqrt_upper(P)+sqrt_upper(barrier)*rootH)**2
    check(budget < old_budget, "strict gain over constant-barrier quadratic: "+label)
    return {
        "label": label, "amplitude": str(amplitude), "chi_radius": str(radius),
        "critical_barrier": str(barrier), "delta": str(delta),
        "reference_budget_upper": str(P), "residual_integral_upper": str(R),
        "F_barrier_upper": str(fb), "F_initial_upper": str(fr),
        "forcing_primitive_cost_upper": str(C),
        "error_X_norm_upper": str(rootQ),
        "solution_X_squared_upper": str(budget),
        "solution_X_squared_for_orientation": float(budget),
        "integer_NSB_exponent": exponent,
        "constant_barrier_budget_for_orientation": float(old_budget)
    }

rows = []
for row in data["finite_candidate_classes"]:
    if not row["passed"]:
        continue
    amplitude = F(row["amplitude"])
    radius = F(row["chi"])
    B = F(row["selected_barrier"])
    check(F(row["endpoint_upper"]) < F(row["Psi_barrier_lower"]),
          "source endpoint margin: "+str(amplitude)+","+str(radius))
    rows.append(compute_budget("new endpoint class", amplitude, radius, B,
        F(row["P2_upper"]), F(row["residual_total_upper"])))
check(len(rows) == 5, "five passing enlarged endpoint classes")

old = compute_budget("both previously proved classes", "11/10 or 1", F(1,20),
                     F(1), F(213,25), F(3,20))
check(F(old["solution_X_squared_upper"]) < 24,
      "common old-class NSB exponent improves from 81 to 24")

result = {
    "status": "PASS_RETAINED_ENDPOINT_HISTORY_NSB_BUDGETS",
    "source_json": str(SOURCE), "source_sha256": EXPECTED_SHA,
    "scope": "Budget consequences of independently certified endpoint gates; smooth full-support initial perturbations allowed.",
    "formula": "F(s)=-2s-4log(1-s/2); C=F(r)+(B/(1-B/2))*Rtot; G=(sqrt(P)+sqrt(F(B))+sqrt(F(B)+C))^2.",
    "new_classes": rows, "old_classes_common_corollary": old,
    "all_cutoff_bound": "sup_K integral ||grad u_K||_2^4 <= (E0*X0/2)*exp(integer_NSB_exponent)",
    "lower_cutoff_argument": "K<sqrt6: sup X_K<=6 X0 and integral X_K^2<=3 E0 X0; e^L>6 for every listed L.",
    "assertion_count": len(checks), "checks": checks,
    "not_claimed": ["Arbitrary-data regularity", "Picard convergence for enlarged endpoint classes", "Optimal amplitude or critical radius"]
}
Path("enlarged_budget_static_result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"], "assertion_count":len(checks),
    "new_classes":[{k:r[k] for k in ["amplitude","chi_radius","critical_barrier",
                    "solution_X_squared_for_orientation","integer_NSB_exponent",
                    "constant_barrier_budget_for_orientation"]} for r in rows],
    "old_common_budget":old["solution_X_squared_for_orientation"],
    "old_common_exponent":old["integer_NSB_exponent"]},indent=2))
