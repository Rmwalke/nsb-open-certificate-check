"""Fresh finite rational rate certificate from the current passing block data.

Reads a JSON certificate as data. No earlier research script is imported
or executed, and no evolution equation or time quadrature is computed.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json


def main():
    source = Path(__file__).with_name("refined_W_block_gate_result.json")
    raw = source.read_bytes()
    data = json.loads(raw)
    passing = [test for test in data["tests"] if test["passed"]]
    assert len(passing) == 1
    test = passing[0]
    assert test["reference"] == "w3"
    matrix = [[Q(value) for value in row]
              for row in test["contraction_matrix"]]
    radii = [Q(value) for value in test["radii_upper"]]
    q = Q(3, 4)
    checks = []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    check("seventeen components", len(matrix) == len(radii) == 17)
    weights_inverse = []
    for i, row in enumerate(matrix):
        check("triangular nonnegative row " + str(i + 1),
              len(row) == i + 1 and all(value >= 0 for value in row))
        check("diagonal below rate " + str(i + 1), row[i] < q)
        prior = sum((row[j] * weights_inverse[j] for j in range(i)), Q(0))
        d = (1 + prior) / (q - row[i])
        weights_inverse.append(d)
        check("positive resolvent vector " + str(i + 1), d > 0)
        check("resolvent identity " + str(i + 1),
              sum((row[j] * weights_inverse[j] for j in range(i + 1)),
                  Q(0)) == q * d - 1)

    scale = max(r / d for r, d in zip(radii, weights_inverse))
    check("positive initial scale", scale > 0)
    for i, (r, d) in enumerate(zip(radii, weights_inverse)):
        check("initial radius dominated " + str(i + 1), r <= scale * d)
    constant_squared = scale**2 * sum((d*d for d in weights_inverse), Q(0))
    check("global X rate constant below fourteen", constant_squared < 14**2)
    check("strict geometric rate", 0 < q < 1)

    result = {
        "status": "PASS_UNIT_W_PICARD_EXPLICIT_GEOMETRIC_RATE",
        "source_json": source.name,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "scope": "Unit-amplitude W, zero perturbation, ordinary global Picard iteration.",
        "sequence": "u^(0)=w2, u^(1)=w3, u^(n+1)=S W-Q(u^(n),u^(n))",
        "rational_rate": str(q),
        "inverse_weights": [str(d) for d in weights_inverse],
        "initial_scale": str(scale),
        "exact_global_constant_squared": str(constant_squared),
        "constant_squared_decimal_orientation_only": float(constant_squared),
        "certified_global_constant": 14,
        "proof": "A d=q d-1<=q d, rho<=scale*d; hence A^j rho<=scale*q^j*d.",
        "rate_from_w3": "||T_W^j(w3)-u||_X <=14*(3/4)^j for every integer j>=0",
        "rate_from_w2": "||u^(n)-u||_X <=14*(3/4)^(n-1) for every integer n>=1",
        "assertion_count": len(checks),
        "checks": checks,
        "not_claimed": [
            "This does not certify a perturbation neighborhood at the same rate.",
            "This is a certificate for the unprojected global iteration, not a sampled trajectory."
        ]
    }
    output = Path(__file__).with_name("W_picard_rate_static_result.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "assertions": len(checks),
        "rate": str(q),
        "global_constant": 14,
        "constant_squared_orientation": float(constant_squared)
    }, indent=2))


if __name__ == "__main__":
    main()
