"""Fresh exact finite-matrix certificates for perturbed heat-seeded Picard.

Only JSON data are read. No previous program is imported or executed.
There is no time evolution or numerical quadrature.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json


def main():
    root = Path(__file__).parent
    names = ["perturbation_block_gate_result.json",
             "perturbed_picard_prefix_static_result.json"]
    raw = {name: (root / name).read_bytes() for name in names}
    hashes = {name: hashlib.sha256(value).hexdigest()
              for name, value in raw.items()}
    gate, prefix = [json.loads(raw[name]) for name in names]
    checks = []

    def check(label, condition):
        assert condition, label
        checks.append(label)

    def integer_sqrt_upper(value):
        check("nonnegative squared constant", value >= 0)
        floor = isqrt(value.numerator // value.denominator)
        return floor if floor * floor >= value else floor + 1

    check("prefix bound to frozen gate",
          prefix["source_sha256"] == hashes[names[0]])
    check("gate passed",
          gate["status"] == "PASS_PERTURBED_W3_BLOCK_CONVERGENCE")
    check("actual heat seed entry passed",
          prefix["status"] ==
          "PASS_PERTURBED_TRUE_HEAT_SEED_FINITE_PREFIX_ENTRY")
    check("common perturbation radius",
          gate["selected_chi"] == prefix["chi_radius"] == "1/10000")
    cert = gate["selected_certificate"]
    n = len(cert["radii_upper"])
    check("seventeen blocks", n == 17)
    c = F(gate["c_upper"])
    V = list(map(F, gate["reference_block_caps"]))
    R = list(map(F, cert["radii_upper"]))
    K = [list(map(F, row)) for row in gate["source_memory_factors"]]
    A = [list(map(F, row)) for row in cert["contraction_matrix"]]
    source = list(map(F, gate["heat_profile_caps"]))
    E3 = list(map(F, prefix["prefix_difference_caps"][3]))
    E = list(map(F, prefix["initial_fixed_point_error_caps"]))
    for i in range(n):
        check(f"matrix row structure {i+1}",
              len(A[i]) == i+1 and all(x >= 0 for x in A[i]))
        check(f"actual seed inside rectangle {i+1}", 0 <= E3[i] < R[i])
        check(f"complete initial error {i+1}", E[i] == E3[i]+R[i])
        check(f"positive source coefficient {i+1}", source[i] > 0)
        for j in range(i+1):
            check(f"matrix formula {i+1},{j+1}",
                  A[i][j] == 2*c*K[i][j]*(V[j]+R[j]))

    # For fixed q, recursively construct the least vector p>=E with Ap<=qp.
    # Positivity and triangularity make each component an explicit maximum.
    rates = []
    for q in (F(3,4), F(4,5), F(9,10)):
        p = []
        for i in range(n):
            check(f"q={q}: diagonal gap {i+1}", A[i][i] < q)
            coupling = sum((A[i][j]*p[j] for j in range(i)), F(0))
            p.append(max(E[i], coupling/(q-A[i][i])))
            check(f"q={q}: initial error dominated {i+1}", p[i] >= E[i] > 0)
            check(f"q={q}: contraction supersolution {i+1}",
                  sum((A[i][j]*p[j] for j in range(i+1)), F(0)) <= q*p[i])
        constant2 = sum((x*x for x in p), F(0))
        cap = integer_sqrt_upper(constant2)
        check(f"q={q}: integer prefactor certified", constant2 <= cap*cap)
        rates.append({
            "q": str(q), "positive_majorant": list(map(str, p)),
            "exact_constant_squared": str(constant2),
            "constant_decimal_orientation_only": float(constant2)**0.5,
            "integer_prefactor": cap,
            "statement": f"||v[n]-u_h||_X <= {cap}*({q})^(n-3), n>=3",
        })

    # Uniform solution-map stability in the same common invariant rectangle.
    # z <= source*chi(h-g)+A z implies z <= b*chi(h-g), by forward substitution.
    b = []
    for i in range(n):
        check(f"stability diagonal gap {i+1}", A[i][i] < 1)
        b.append((source[i]+sum((A[i][j]*b[j] for j in range(i)), F(0)))
                 /(1-A[i][i]))
        check(f"positive stability coefficient {i+1}", b[i] > 0)
        check(f"stability inverse identity {i+1}",
              b[i]-sum((A[i][j]*b[j] for j in range(i+1)), F(0)) == source[i])
    lipschitz2 = sum((x*x for x in b), F(0))
    lipschitz = integer_sqrt_upper(lipschitz2)
    check("integer solution-map Lipschitz bound",
          lipschitz2 <= lipschitz*lipschitz)

    budget = sum(((V[i]+R[i])**2 for i in range(n)), F(0))
    check("stored budget matches", budget == F(cert["solution_X_squared_upper"]))
    check("uniform budget below nine", budget < 9)
    # e^9 >= 1+9+9^2/2 >24 suffices to absorb every cutoff K<sqrt(24).
    check("small-cutoff absorption", 1+F(9)+F(81,2) > 24)

    result = {
        "status": "PASS_PERTURBED_PICARD_RATES_AND_SOLUTION_STABILITY",
        "source_sha256": hashes,
        "scope": "Smooth real mean-zero divergence-free h,g with A_-1(h),A_-1(g)<=1/10000.",
        "sequence": "v[0]=S(W+h), v[n+1]=S(W+h)-Q(v[n],v[n])",
        "rates": rates,
        "rate_proof": "E=E3+R<=p and Ap<=q*p imply A^j E<=q^j*p, hence the l2 prefactor.",
        "solution_stability": {
            "source_coefficients": list(map(str, source)),
            "inverse_vector": list(map(str, b)),
            "exact_Lipschitz_constant_squared": str(lipschitz2),
            "Lipschitz_decimal_orientation_only": float(lipschitz2)**0.5,
            "integer_Lipschitz_constant": lipschitz,
            "statement": f"||u_h-u_g||_X <= {lipschitz}*A_-1(h-g)",
            "proof": "Both solutions lie in the common rectangle. Their difference satisfies z<=source*A_-1(h-g)+Az; forward substitution gives z<=(I-A)^(-1)source*A_-1(h-g).",
        },
        "solution_X_squared_upper": str(budget),
        "all_cutoff_NSB": "(E0*X0/2)*exp(9), E0=||W+h||_2^2, X0=||grad(W+h)||_2^2",
        "large_cutoff_budget": "K>=sqrt(24): integral A0(u_K)^2 <9",
        "small_cutoff_argument": "K<sqrt(24): integral ||grad u_K||_2^4 <=K^2*E0^2/2 <=12*E0*X0 <(E0*X0/2)*exp(9).",
        "assertion_count": len(checks),
        "checks": checks,
        "limits": [
            "Rates hold in X, not asserted in H1.",
            "The neighborhood is contained in an earlier larger regularity class.",
            "Failure of a tested larger radius does not show iteration divergence.",
            "No arbitrary-data theorem or external/formal verification claimed.",
        ],
    }
    destination = root / "perturbed_stability_rate_static_result.json"
    destination.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "assertions": len(checks),
        "rates": [{"q": item["q"], "prefactor": item["integer_prefactor"],
                   "constant_orientation": item["constant_decimal_orientation_only"]}
                  for item in rates],
        "solution_Lipschitz_constant": lipschitz,
        "Lipschitz_orientation": float(lipschitz2)**0.5,
        "budget_orientation": float(budget),
    }, indent=2))


if __name__ == "__main__":
    main()
