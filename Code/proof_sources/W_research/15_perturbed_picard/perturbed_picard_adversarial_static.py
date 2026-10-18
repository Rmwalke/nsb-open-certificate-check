"""Fresh data-only audit of the perturbation block certificate.

No imports or execution of any earlier research source; no dynamics or
quadrature. The frozen unit-W certificate is assumed as already audited.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import hashlib
import json


@lru_cache(None)
def exponential_interval(x):
    assert x >= 0
    n = max(100, 4*(x.numerator//x.denominator+1))
    terms = [x**j/factorial(j) for j in range(n+1)]
    low = sum(terms, F(0))
    high = low + terms[-1]*x/(n+1)/(1-x/(n+2))
    return 1/high, 1/low


def heat_factor_squared_upper(left, length):
    if left == 0:
        return F(1)
    candidates = [exponential_interval(2*left)[1]]
    if length is not None:
        x = length/(2*left+length)
        n = 48
        loglo = 2*sum((x**(2*j+1)/(2*j+1) for j in range(n)), F(0))
        loghi = loglo + 2*x**(2*n+1)/((2*n+1)*(1-x*x))
        candidates.append(length/(left+length)*exponential_interval(left/length*loglo)[1])
        if loghi <= 2*length:
            candidates.append(exponential_interval(2*left)[1]*(1-exponential_interval(2*length)[0]))
    return min(candidates)


def main():
    names = ['refined_W_block_gate_result.json', 'perturbation_block_gate_result.json',
             'perturbed_picard_prefix_static_result.json',
             'perturbed_stability_rate_static_result.json']
    blobs = {name: Path(name).read_bytes() for name in names}
    hashes = {name: hashlib.sha256(raw).hexdigest() for name, raw in blobs.items()}
    base, current, prefix, rates = [json.loads(raw) for raw in blobs.values()]
    assert hashes[names[0]] == '17b578bd37cb971cfbd21a58749bc95a94d62e126fcced0fe21bb3b85a3420b0'
    assert hashes[names[1]] == 'c789c350a041c8e9eba01ccea50f1839fbe2a678439f924f6c13e4f29f2ca461'
    assert current['source_sha256'] == hashes[names[0]]
    source = next(t for t in base['tests'] if t['passed'])
    cert = current['selected_certificate']
    epsilon = F(current['selected_chi'])
    assert epsilon == F(1,10000) == F(cert['chi'])
    c = F(current['c_upper'])
    assert c > 0 and c*c > F(1,8)
    V = [F(x) for x in current['reference_block_caps']]
    D = [F(x) for x in current['baseline_defect_block_caps']]
    K = [[F(x) for x in row] for row in current['source_memory_factors']]
    H = [F(x) for x in current['heat_profile_caps']]
    rho = [F(x) for x in cert['radii_upper']]
    A = [[F(x) for x in row] for row in cert['contraction_matrix']]
    assert all(len(a) == 17 for a in [V,D,K,H,rho,A,cert['rows']])
    slacks = []
    for i, row in enumerate(cert['rows']):
        sr = source['rows'][i]
        assert V[i] == F(sr['V_upper']) and D[i] == F(sr['defect_upper'])
        assert K[i] == [F(x) for x in sr['coefficients']]+[F(1)]
        left = F(sr['left'])
        length = None if sr['right'] == 'infinity' else F(sr['right'])-left
        assert H[i] > 0
        # H_i is the coefficient per unit chi. Compare to the exact 1/sqrt(2),
        # rather than relying on the original program's selected c upper.
        assert H[i]*H[i] >= heat_factor_squared_upper(left,length)/2
        assert F(row['initial_heat_cap']) == epsilon*H[i]
        assert F(row['perturbed_defect_cap']) == D[i]+epsilon*H[i]
        assert F(row['radius_upper']) == rho[i] and rho[i] > 0
        forcing = D[i]+epsilon*H[i]
        nonlinear = c*sum((K[i][j]*(2*V[j]*rho[j]+rho[j]**2) for j in range(i+1)), F(0))
        slack = rho[i]-forcing-nonlinear
        assert slack > 0
        slacks.append(slack)
        assert len(A[i]) == i+1
        for j in range(i+1):
            assert A[i][j] == 2*c*K[i][j]*(V[j]+rho[j])
        assert A[i][i] < F(3,5)
    spectral = max(A[i][i] for i in range(17))
    assert F(cert['spectral_radius_upper']) == spectral
    budget = sum(((V[i]+rho[i])**2 for i in range(17)), F(0))
    assert budget == F(cert['solution_X_squared_upper']) < 9
    assert prefix['source_sha256'] == hashes[names[1]]
    assert prefix['baseline_source_sha256'] == hashes[names[0]]
    assert F(prefix['chi_radius']) == epsilon
    assert prefix['grid'] == current['grid'] == source['grid']
    references = [[F(x) for x in row] for row in prefix['unperturbed_reference_caps']]
    differences = [[F(x) for x in row] for row in prefix['prefix_difference_caps']]
    assert len(references) == 3 and len(differences) == 4

    def reference_tail(t, order):
        if t is None:
            return F(0), F(0)
        terms = [(F(9,2),2)]
        if order == 1:
            terms += [(6*t+2,3), (F(9,4)*t*t+F(9,8)*t+F(9,32),4)]
        lo = hi = F(0)
        for coefficient, rate in terms:
            a,b = exponential_interval(rate*t)
            lo += coefficient*a
            hi += coefficient*b
        return lo,hi

    for i,sr in enumerate(source['rows']):
        left = F(sr['left'])
        right = None if sr['right'] == 'infinity' else F(sr['right'])
        for order in (0,1):
            moment_cap = reference_tail(left,order)[1]-reference_tail(right,order)[0]
            assert references[order][i]**2 >= moment_cap
        assert references[2][i] == F(sr['reference2_norm_upper'])
        assert differences[0][i] == epsilon*H[i]
    for order in range(3):
        for i in range(17):
            bound = epsilon*H[i] + c*sum((K[i][j]*(2*references[order][j]*differences[order][j]
                       + differences[order][j]**2) for j in range(i+1)), F(0))
            assert differences[order+1][i] >= bound
    for i in range(17):
        assert 0 <= differences[3][i] < rho[i]
        assert F(prefix['entry_margins'][i]) == rho[i]-differences[3][i]
        assert F(prefix['initial_fixed_point_error_caps'][i]) == rho[i]+differences[3][i]
    max_entry = max(differences[3][i]/rho[i] for i in range(17))
    assert max_entry == F(prefix['maximum_entry_ratio']) < 1
    assert rates['source_sha256'][names[1]] == hashes[names[1]]
    assert rates['source_sha256'][names[2]] == hashes[names[2]]
    E = [rho[i]+differences[3][i] for i in range(17)]
    verified_rates = []
    for item in rates['rates']:
        q = F(item['q'])
        p = [F(x) for x in item['positive_majorant']]
        assert len(p) == 17 and 0 < q < 1
        for i in range(17):
            assert p[i] >= E[i] > 0
            assert sum((A[i][j]*p[j] for j in range(i+1)),F(0)) <= q*p[i]
        constant_squared = sum((x*x for x in p),F(0))
        assert F(item['exact_constant_squared']) == constant_squared <= item['integer_prefactor']**2
        verified_rates.append({'q':str(q),'prefactor':item['integer_prefactor']})
    assert {'q':'9/10','prefactor':6} in verified_rates
    stability = rates['solution_stability']
    b = [F(x) for x in stability['inverse_vector']]
    assert [F(x) for x in stability['source_coefficients']] == H
    for i in range(17):
        assert b[i] > 0
        assert b[i]-sum((A[i][j]*b[j] for j in range(i+1)),F(0)) == H[i]
    stability_squared = sum((x*x for x in b),F(0))
    assert F(stability['exact_Lipschitz_constant_squared']) == stability_squared <= 4757**2
    assert stability['integer_Lipschitz_constant'] == 4757
    assert F(rates['solution_X_squared_upper']) == budget
    result = {
        'status': 'PASS_PERTURBED_TRUE_HEAT_SEED_INDEPENDENT_AUDIT',
        'source_sha256': hashes,
        'chi_radius': str(epsilon),
        'heat_caps_independently_checked': 17,
        'strict_invariance_rows': 17,
        'all_new_matrix_entries_checked': True,
        'all_diagonals_below_3_over_5': True,
        'spectral_upper_decimal_orientation_only': float(spectral),
        'solution_X_squared_upper_below_9': True,
        'solution_X_squared_decimal_orientation_only': float(budget),
        'minimum_positive_invariance_slack_orientation_only': float(min(slacks)),
        'unperturbed_prefix_reference_norm_caps_independently_checked': 51,
        'perturbed_finite_prefix_update_caps_checked': 51,
        'strict_third_iterate_entry_blocks': 17,
        'maximum_entry_ratio_decimal_orientation_only': float(max_entry),
        'finite_prefix_entry': 'Actual perturbed heat-seeded third iterate lies strictly inside the rectangle.',
        'independently_verified_geometric_rates': verified_rates,
        'independently_verified_solution_map_Lipschitz_constant': 4757,
        'scope': 'New iteration convergence inside the previously regular chi<=1/20 class.'
    }
    Path('perturbed_picard_adversarial_static_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
