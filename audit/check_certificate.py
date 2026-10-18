#!/usr/bin/env python3
"""Data-only exact-rational audit of manuscript Section 2.5.

No original proof program is imported or executed. Archived analytic
enclosures are INPUT ASSUMPTIONS, not established by this checker.
Usage: python audit_w15.py --code-dir /path/to/Code --output-dir /path/to/out
"""
import argparse
import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def run(code_dir, out_dir):
    root = Path(code_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        'gate': 'proof_sources/W_research/15_perturbed_picard/perturbation_block_gate_result.json',
        'prefix': 'proof_sources/W_research/15_perturbed_picard/perturbed_picard_prefix_static_result.json',
        'base': 'proof_sources/W_research/14_unit_W_picard/refined_W_block_gate_result.json',
    }
    docs, inputs = {}, {}
    for label, rel in paths.items():
        raw = (root / rel).read_bytes()
        docs[label] = json.loads(raw)
        inputs[label] = {'relative_path': rel, 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
    g, p, b = docs['gate'], docs['prefix'], docs['base']
    checks = []
    def check(name, ok):
        checks.append({'name': name, 'passed': bool(ok)})
        if not ok:
            raise ValueError('Input consistency check failed: ' + name)

    check('gate binds baseline bytes', g['source_sha256'] == inputs['base']['sha256'])
    check('prefix binds gate bytes', p['source_sha256'] == inputs['gate']['sha256'])
    check('prefix binds baseline bytes', p['baseline_source_sha256'] == inputs['base']['sha256'])
    grid = g['grid']
    check('17 intervals including final infinite interval', len(grid) == 18 and grid[0] == '0' and grid[-1] == 'infinity')
    check('finite endpoints increasing', all(F(grid[i]) < F(grid[i+1]) for i in range(16)))
    check('prefix grid agrees', p['grid'] == grid)
    selected = g['selected_certificate']
    n = 17
    c, eps = F(g['c_upper']), F(g['selected_chi'])
    check('radius exactly 1/10000', eps == F(1, 10000) == F(p['chi_radius']) == F(selected['chi']))
    check('positive valid bilinear constant majorant', c > 0 and c*c >= F(1, 8))
    V = list(map(F, g['reference_block_caps']))
    d = list(map(F, g['baseline_defect_block_caps']))
    s = list(map(F, g['heat_profile_caps']))
    rho = list(map(F, selected['radii_upper']))
    K = [list(map(F, row)) for row in g['source_memory_factors']]
    Vpre = [list(map(F, row)) for row in p['unperturbed_reference_caps']]
    Eold = [list(map(F, row)) for row in p['prefix_difference_caps']]
    check('vector dimensions', all(len(v) == n for v in [V, d, s, rho, *Vpre, *Eold]) and len(Vpre) == 3 and len(Eold) == 4 and len(K) == n)
    check('nonnegative input bounds, positive radii', all(x >= 0 for v in [V,d,s,*Vpre,*Eold] for x in v) and all(x > 0 for x in rho))
    base_rows = [x for x in b['tests'] if x['name'] == 'w3_prior17_exact_U2_plus_D3']
    check('unique baseline block collection', len(base_rows) == 1 and len(base_rows[0]['rows']) == n)
    check('baseline grid agrees', base_rows[0]['grid'] == grid)
    for i, row in enumerate(base_rows[0]['rows']):
        check(f'block {i+1}: triangular complete kernel row', len(K[i]) == i+1 and K[i][-1] == 1 and all(x >= 0 for x in K[i]))
        check(f'block {i+1}: prior kernel factors bound to baseline', K[i][:-1] == list(map(F, row['coefficients'])))
        check(f'block {i+1}: reference cap bound to baseline', V[i] == F(row['V_upper']))
        check(f'block {i+1}: defect cap bound to baseline', d[i] == F(row['defect_upper']))
        check(f'block {i+1}: second-iterate cap bound to baseline', Vpre[2][i] == F(row['reference2_norm_upper']))
        check(f'block {i+1}: initial perturbation bookkeeping', eps*s[i] == F(selected['rows'][i]['initial_heat_cap']) == Eold[0][i])

    # Reconstruct the three-step majorant with exact fractions and no rounding.
    Enew = [[eps*x for x in s]]
    recurrence_margins = []
    for stage in range(3):
        new, archived_margins = [], []
        for i in range(n):
            def rhs(values):
                return eps*s[i] + c*sum((K[i][j]*(2*Vpre[stage][j]*values[j] + values[j]**2) for j in range(i+1)), F(0))
            new.append(rhs(Enew[stage]))
            archived_margins.append(Eold[stage+1][i] - rhs(Eold[stage]))
        Enew.append(new)
        recurrence_margins.append(archived_margins)

    rows = []
    for i in range(n):
        lhs = d[i] + eps*s[i] + c*sum((K[i][j]*(2*V[j]*rho[j]+rho[j]**2) for j in range(i+1)), F(0))
        inv_margin = rho[i] - lhs
        entry_margin = rho[i] - Enew[3][i]
        old_entry_margin = rho[i] - Eold[3][i]
        rec_ok = all(stage[i] >= 0 for stage in recurrence_margins)
        rows.append({
            'block': i+1, 'left': grid[i], 'right': grid[i+1],
            'radius_exact': str(rho[i]), 'invariance_lhs_exact': str(lhs),
            'invariance_margin_exact': str(inv_margin),
            'invariance_margin_display': format(float(inv_margin), '.12e'),
            'invariance_pass': inv_margin >= 0,
            'fresh_E3_exact': str(Enew[3][i]), 'fresh_entry_margin_exact': str(entry_margin),
            'fresh_entry_ratio_display': format(float(Enew[3][i]/rho[i]), '.15g'),
            'fresh_entry_pass': entry_margin > 0,
            'archived_E3_exact': str(Eold[3][i]), 'archived_entry_margin_exact': str(old_entry_margin),
            'archived_entry_pass': old_entry_margin > 0,
            'archived_recurrences_pass': rec_ok,
            'fresh_le_archived': Enew[3][i] <= Eold[3][i],
        })
    passed = all(all(row[k] for k in ['invariance_pass','fresh_entry_pass','archived_entry_pass','archived_recurrences_pass','fresh_le_archived']) for row in rows)
    result = {
        'status': 'PASS_DATA_ONLY' if passed else 'FAIL_DATA_ONLY',
        'scope': 'Exact rational arithmetic implications of supplied archived enclosures for 17-block invariance and third-iterate entry only.',
        'assumptions_and_limits': [
            'Input V, d, s, K and unperturbed prefix caps are assumed to be mathematically valid enclosures.',
            'Hashes and equality checks establish consistency of these files, not validity or provenance of their mathematical claims.',
            'No original scientific script, solver or verifier was run or imported.',
            'No independent Fourier reconstruction, kernel enclosure proof, PDE limit argument, novelty assessment or peer review is supplied.',
            'Displayed decimals are informational; every pass/fail decision uses exact Fraction comparisons.',
            'The three-step majorant is recomputed without rounding; archived prefix recurrence inequalities are separately checked.',
        ],
        'inputs': inputs, 'input_consistency_checks': checks,
        'formula_reference': 'manuscript/sections/02_main_proof.tex: main:invariance, main:prefix, main:accepted',
        'fresh_prefix_caps_exact': [[str(x) for x in v] for v in Enew],
        'archived_recurrence_margins_exact': [[str(x) for x in v] for v in recurrence_margins],
        'summary': {
            'blocks': n, 'invariance_pass_count': sum(r['invariance_pass'] for r in rows),
            'entry_pass_count': sum(r['fresh_entry_pass'] for r in rows),
            'archived_prefix_recurrence_pass_count': sum(x >= 0 for stage in recurrence_margins for x in stage),
            'smallest_invariance_margin_exact': str(min(F(r['invariance_margin_exact']) for r in rows)),
            'smallest_invariance_margin_display': format(float(min(F(r['invariance_margin_exact']) for r in rows)), '.16e'),
            'largest_fresh_entry_ratio_display': format(float(max(Enew[3][i]/rho[i] for i in range(n))), '.16g'),
        }, 'blocks': rows,
    }
    (out/'W15_data_only_audit.json').write_text(json.dumps(result, indent=2)+'\n')
    with (out/'W15_blocks.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({'status': result['status'], **result['summary']}, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--code-dir', required=True)
    ap.add_argument('--output-dir', required=True)
    args = ap.parse_args()
    raise SystemExit(run(args.code_dir, args.output_dir))
