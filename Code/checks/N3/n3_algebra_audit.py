"""Independent exact-rational Fourier audit; no time integration.

Reads only n3_exact_algebra_intermediate.json and writes only
n3_algebra_audit_result.json in this directory.
"""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path

BASE = Path(__file__).parent
SOURCE = BASE / 'n3_exact_algebra_intermediate.json'
OUTPUT = BASE / 'n3_algebra_audit_result.json'
Z = Q(0)


def clean(field):
    return {key: value for key, value in field.items() if value}


def read_field(entries):
    field = defaultdict(Q)
    seen = set()
    duplicates = zeros = 0
    for entry in entries:
        key = (tuple(entry['k']), entry['c'])
        duplicates += key in seen
        seen.add(key)
        value = Q(entry['v'])
        zeros += not value
        field[key] += value
    return clean(field), duplicates, zeros


def vectors(field):
    out = {}
    for (k, component), value in field.items():
        out.setdefault(k, [Z, Z, Z])[component] = value
    return out


def inner(left, right):
    return sum((value * right.get(key, Z) for key, value in left.items()), Z)


def add(left, right):
    out = defaultdict(Q, left)
    for key, value in right.items():
        out[key] += value
    return clean(out)


def mismatch(left, right):
    return sum(left.get(key, Z) != right.get(key, Z)
               for key in left.keys() | right.keys())


def shell(k):
    return sum(component * component for component in k)


def raw_bilinear(left, right):
    # Fourier coefficients are i*a and i*b. Thus raw B/i = -(a_p dot q)*b_q.
    out = defaultdict(Q)
    for p, a in vectors(left).items():
        for q, b in vectors(right).items():
            factor = -sum((a[c] * q[c] for c in range(3)), Z)
            if factor:
                k = tuple(p[c] + q[c] for c in range(3))
                for c in range(3):
                    if b[c]:
                        out[(k, c)] += factor * b[c]
    return clean(out)


def leray(raw):
    out, zero_mode = {}, {}
    for k, vector in vectors(raw).items():
        radius2 = shell(k)
        if not radius2:
            for c in range(3):
                if vector[c]:
                    zero_mode[(k, c)] = vector[c]
            continue
        dot = sum((k[c] * vector[c] for c in range(3)), Z)
        for c in range(3):
            value = vector[c] - Q(k[c], radius2) * dot
            if value:
                out[(k, c)] = value
    return out, zero_mode


def field_checks(field):
    return {
        'divergence': sum(sum((k[c] * v[c] for c in range(3)), Z) != 0
                          for k, v in vectors(field).items()),
        'conjugacy': sum(field.get((tuple(-z for z in k), c), Z) != -v
                         for (k, c), v in field.items()),
        'zero_mode': sum(k == (0, 0, 0) for k, c in field),
    }


def audit():
    source_bytes = SOURCE.read_bytes()
    data = json.loads(source_bytes)
    basis, checks = [], Counter()
    for entries in data['basis']:
        field, duplicates, zeros = read_field(entries)
        basis.append(field)
        checks['duplicate_field_coordinates'] += duplicates
        checks['explicit_zero_field_coefficients'] += zeros
    count = len(basis)
    checks['basis_count_not_ten'] = count != 10
    norms = [inner(field, field) for field in basis]
    for i, field in enumerate(basis):
        checks.update({'basis_' + key: value for key, value in field_checks(field).items()})
        checks['shell_assignment'] += sum(shell(k) != data['shells'][i] for k, c in field)
        checks['norm'] += norms[i] != Q(data['norms'][i])
        checks['zero_basis_norm'] += norms[i] == 0
    checks['offdiagonal_inner_products'] = sum(
        inner(basis[i], basis[j]) != 0
        for i in range(count) for j in range(i + 1, count))

    stored_q = {}
    for output, i, j, value in data['ode_coefficients']:
        key = (output, i, j)
        checks['duplicate_ode_keys'] += key in stored_q
        stored_q[key] = Q(value)
    stored_r = {}
    for item in data['residual_pairs']:
        key = (item['i'], item['j'])
        checks['duplicate_residual_pair_keys'] += key in stored_r
        field, duplicates, zeros = read_field(item['fourier'])
        checks['duplicate_field_coordinates'] += duplicates
        checks['explicit_zero_field_coefficients'] += zeros
        stored_r[key] = field

    valid_pairs = list(combinations_with_replacement(range(count), 2))
    valid_q_keys = {(out, i, j) for out in range(count) for i, j in valid_pairs}
    checks['unexpected_ode_keys'] = len(set(stored_q) - valid_q_keys)
    checks['unexpected_residual_pair_keys'] = len(set(stored_r) - set(valid_pairs))
    derived_q, exterior_pairs, failed_pairs = {}, {}, []
    for i, j in valid_pairs:
        raw = raw_bilinear(basis[i], basis[j])
        if i != j:
            raw = add(raw, raw_bilinear(basis[j], basis[i]))
        projected, zero_mode = leray(raw)
        checks['raw_zero_mode_coefficients'] += len(zero_mode)
        checks.update({'interaction_' + key: value
                       for key, value in field_checks(projected).items()})
        retained = {key: value for key, value in projected.items() if shell(key[0]) <= 9}
        exterior = {key: value for key, value in projected.items() if shell(key[0]) > 9}
        derived_reconstruction, stored_reconstruction = defaultdict(Q), defaultdict(Q)
        q_errors = 0
        for out in range(count):
            coefficient = inner(retained, basis[out]) / norms[out]
            if coefficient:
                derived_q[(out, i, j)] = coefficient
            stored_coefficient = stored_q.get((out, i, j), Z)
            q_errors += coefficient != stored_coefficient
            for key, value in basis[out].items():
                derived_reconstruction[key] += coefficient * value
                stored_reconstruction[key] += stored_coefficient * value
        closure_errors = mismatch(retained, clean(derived_reconstruction))
        retained_errors = mismatch(retained, clean(stored_reconstruction))
        exterior_errors = mismatch(exterior, stored_r.get((i, j), {}))
        checks['ode_coefficients'] += q_errors
        checks['closure_coefficients'] += closure_errors
        checks['closure_pairs'] += bool(closure_errors)
        checks['retained_reconstruction_coefficients'] += retained_errors
        checks['exterior_coefficients'] += exterior_errors
        if q_errors or closure_errors or retained_errors or exterior_errors:
            failed_pairs.append([i, j, q_errors, closure_errors, retained_errors, exterior_errors])
        if exterior:
            exterior_pairs[(i, j)] = exterior

    energy = defaultdict(Q)
    for (out, i, j), coefficient in derived_q.items():
        energy[tuple(sorted((out, i, j)))] += norms[out] * coefficient
    checks['energy_cubic_coefficients'] = len(clean(energy))

    w, d = {}, {}
    for c, axis in [(0, 2), (1, 0), (2, 1)]:
        for sign in [-1, 1]:
            k = [0, 0, 0]
            k[axis] = sign
            w[(tuple(k), c)] = Q(-sign, 2)
    for c, sine_axis, double_axis in [(0, 2, 1), (1, 0, 2), (2, 1, 0)]:
        for a in [-1, 1]:
            for b in [-1, 1]:
                k = [0, 0, 0]
                k[sine_axis], k[double_axis] = a, 2 * b
                d[(tuple(k), c)] = Q(-a, 8)
    c, unused = leray(raw_bilinear(w, w))
    f, unused = leray(raw_bilinear(c, c))
    for name, expected, index in [('W', w, 0), ('C', c, 1), ('D', d, 2), ('F', f, 4)]:
        checks['named_field_' + name] = mismatch(expected, basis[index])

    result = {
        'status': 'PASS_EXACT_ALGEBRA' if not any(checks.values()) else 'FAIL_EXACT_ALGEBRA',
        'source_sha256': hashlib.sha256(source_bytes).hexdigest(),
        'arithmetic': 'Python Fraction, exact; independent Fourier convolution implementation',
        'fourier_convention': 'u_k=i*a_k; stored quadratic coefficients are +P3B; ODE subtracts them',
        'mismatch_counts': dict(sorted(checks.items())),
        'basis_count': count, 'normalized_squared_norms': [str(value) for value in norms],
        'shell_dimensions': dict(sorted(Counter(data['shells']).items())),
        'retained_signed_modes': len({k for field in basis for k, c in field}),
        'pairs_checked': len(valid_pairs), 'nonzero_ode_coefficients': len(derived_q),
        'nonzero_exterior_pairs': len(exterior_pairs),
        'exterior_signed_modes': len({k for field in exterior_pairs.values() for k, c in field}),
        'exterior_squared_shells': sorted({shell(k) for field in exterior_pairs.values() for k, c in field}),
        'WD_exterior_squared_norm': str(inner(exterior_pairs[(0, 2)], exterior_pairs[(0, 2)])),
        'failed_pairs': failed_pairs,
    }
    OUTPUT.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if any(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    audit()
