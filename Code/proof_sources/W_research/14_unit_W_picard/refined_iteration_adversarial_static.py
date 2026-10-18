"""Independent rational audit of supplied w3 block certificate, data only.

Does not import or execute any prior research program. No evolution or
quadrature. This checks the saved certificate and its source bindings.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from functools import lru_cache
import hashlib
import json


@lru_cache(None)
def exp_negative_interval(x):
    """Positive Taylor sum for exp(x), then reciprocal; exact rationals."""
    assert x >= 0
    n = max(100, 4 * (x.numerator // x.denominator + 1))
    terms = [x**j / factorial(j) for j in range(n + 1)]
    lower = sum(terms, F(0))
    next_term = terms[-1] * x / (n + 1)
    upper = lower + next_term / (1 - x / (n + 2))
    return 1 / upper, 1 / lower


def valid_kernel_upper(kappa, gap, length):
    if gap == 0:
        return kappa >= 1
    bounds = [exp_negative_interval(2 * gap)[1]]
    if length is not None:
        x = length / (2 * gap + length)
        n = 48
        log_lower = 2 * sum((x**(2*j+1) / (2*j+1) for j in range(n)), F(0))
        log_upper = log_lower + 2*x**(2*n+1) / ((2*n+1)*(1-x*x))
        bounds.append(length / (gap+length) *
                      exp_negative_interval(gap / length * log_lower)[1])
        if log_upper <= 2 * length:
            bounds.append(exp_negative_interval(2*gap)[1] *
                          (1-exp_negative_interval(2*length)[0]))
    return kappa*kappa >= min(bounds)


def main():
    files = {name: Path(name).read_bytes() for name in (
        'third_heat_defect_static_result.json',
        'fourth_heat_defect_static_result.json',
        'refined_W_block_gate_result.json')}
    third, fourth, full = [json.loads(v) for v in files.values()]
    hashes = {name: hashlib.sha256(v).hexdigest() for name, v in files.items()}
    assert full['source_sha256'] == hashes['third_heat_defect_static_result.json']
    assert fourth['source']['sha256'] == hashes['third_heat_defect_static_result.json']
    assert full['fourth_source_sha256'] == hashes['fourth_heat_defect_static_result.json']
    assert fourth['definitions']['amplitude_a'] == '1/1'
    assert max(F(m['lambda']) for m in third['all_defect_modes']) <= 22
    decrease = fourth['actual_single_step_decrease']
    l2_squared = sum((F(m['exact_L2_time_squared']) for m in third['all_defect_modes']), F(0))
    assert l2_squared == F(decrease['z3_spacetime_L2_exact_squared'])
    assert F(decrease['z3_X_rational_lower_bound'])**2 <= l2_squared
    assert F(decrease['z4_X_rational_upper_bound']) < F(decrease['z3_X_rational_lower_bound'])
    tests = [t for t in full['tests'] if t['passed']]
    assert len(tests) == 1
    test = tests[0]
    assert test['reference'] == 'w3' and not test['refined_diagonal']
    rows = test['rows']
    assert len(rows) == len(fourth['block_defect_bounds']) == 17
    c = F(707106781187, 2000000000000)
    assert c*c > F(1,8)
    radii = [F(row['radius_upper']) for row in rows]
    V = [F(row['V_upper']) for row in rows]
    A = [[F(x) for x in row] for row in test['contraction_matrix']]
    slacks, diagonal_slacks = [], []
    kernel_checks = 0
    for i, row in enumerate(rows):
        assert row['status'] == 'PASS'
        source = fourth['block_defect_bounds'][i]
        assert F(row['left']) == F(source['left'])
        if row['right'] == 'infinity':
            assert source['right'] == 'infinity' and i == 16
            length = None
        else:
            assert F(row['right']) == F(source['right'])
            length = F(row['right']) - F(row['left'])
        defect = F(row['defect_upper'])
        assert defect == F(source['full_defect_norm_upper'])
        assert defect >= F(source['mixed_defect_norm_upper'])
        assert V[i] == F(row['reference2_norm_upper']) + F(row['third_heat_defect_upper'])
        assert radii[i] > 0
        kernels = [F(x) for x in row['coefficients']]
        assert len(kernels) == i and len(A[i]) == i+1
        history = F(0)
        for j, k in enumerate(kernels):
            gap = F(row['left']) - F(rows[j]['right'])
            assert valid_kernel_upper(k, gap, length)
            kernel_checks += 1
            history += c*k*(2*V[j]*radii[j] + radii[j]**2)
            assert A[i][j] == 2*c*k*(V[j]+radii[j])
        assert history <= F(row['history_upper'])
        diagonal = 2*c*V[i]
        assert diagonal == F(row['diagonal_linear_upper'])
        assert A[i][i] == diagonal+2*c*radii[i]
        slack = radii[i]-defect-history-diagonal*radii[i]-c*radii[i]**2
        assert slack > 0
        slacks.append(slack)
        diagonal_slacks.append(1-A[i][i])
        assert A[i][i] < F(3,5)
    spectral_upper = max(A[i][i] for i in range(17))
    assert spectral_upper == F(test['spectral_radius_upper'])
    budget = sum(((V[i]+radii[i])**2 for i in range(17)), F(0))
    assert budget < 8
    result = {
        'status': 'PASS_INDEPENDENT_W3_RECTANGLE_AUDIT',
        'source_sha256': hashes,
        'block_count': 17,
        'independently_checked_offdiagonal_kernels': kernel_checks,
        'all_invariance_slacks_strictly_positive': True,
        'all_diagonal_lipschitz_entries_below_3_over_5': True,
        'spectral_radius_upper_decimal_for_orientation': float(spectral_upper),
        'solution_X_squared_upper_below_8': True,
        'solution_X_squared_decimal_for_orientation': float(budget),
        'actual_single_step_decrease_lower_upper_comparison': True,
        'minimum_invariance_slack_positive_decimal_for_orientation': float(min(slacks)),
        'scope': 'Unit W only; certificate data read without executing old code. '
                 'All-past z4 and pointwise quadratic envelope are audited analytically '
                 'in refined_iteration_adversarial.md. Rate calculation is owned by root.'
    }
    Path('refined_iteration_adversarial_static_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
