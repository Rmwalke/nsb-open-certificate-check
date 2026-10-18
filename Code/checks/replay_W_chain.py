"""Replay the finite exact-arithmetic W Stage12--16 certificate chain.

This is an orchestration/receipt script, not a new Fourier checker. It copies
15 explicitly selected archived programs into an isolated directory, generates
their JSON inputs in dependency order, and compares outputs with the immutable
archive. No PDE/ODE time integration or numerical quadrature is requested.

From the code-package root:
  python3 checks/replay_W_chain.py
The default source is resolved relative to this script. The default output is
replay_new_W in the caller's current working directory. Receipt locations use
logical source/output roots and replayable commands, never absolute local paths.
The output directory must be absent or empty. The source is never written.
"""
import argparse
import ast
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

CHECKERS = [
 ('12_second_correction','second_correction_fourier_check.py','second_correction_fourier_check_result.json'),
 ('12_second_correction','second_corrected_gate_static.py','second_corrected_gate_static_result.json'),
 ('13_iteration_convergence','third_heat_defect_static.py','third_heat_defect_static_result.json'),
 ('14_unit_W_picard','fourth_heat_defect_static.py','fourth_heat_defect_static_result.json'),
 ('14_unit_W_picard','refined_W_block_gate_check.py','refined_W_block_gate_result.json'),
 ('14_unit_W_picard','refined_iteration_adversarial_static.py','refined_iteration_adversarial_static_result.json'),
 ('14_unit_W_picard','W_picard_rate_static.py','W_picard_rate_static_result.json'),
 ('15_perturbed_picard','perturbation_block_gate.py','perturbation_block_gate_result.json'),
 ('15_perturbed_picard','perturbed_picard_prefix_static.py','perturbed_picard_prefix_static_result.json'),
 ('15_perturbed_picard','perturbed_stability_rate_static.py','perturbed_stability_rate_static_result.json'),
 ('15_perturbed_picard','perturbed_picard_adversarial_static.py','perturbed_picard_adversarial_static_result.json'),
 ('16_enlarged_endpoint','enlarged_endpoint_gate.py','enlarged_endpoint_gate_result.json'),
 ('16_enlarged_endpoint','enlarged_budget_static.py','enlarged_budget_static_result.json'),
 ('16_enlarged_endpoint','enlarged_endpoint_adversarial_static.py','enlarged_endpoint_adversarial_static_result.json'),
 ('16_enlarged_endpoint','enlarged_class_inclusion_static.py','enlarged_class_inclusion_static_result.json'),
]
CAP_MAP={'F':'W,W','K_WC':'W,G5','K_WH':'W,H','G':'W,F','K_FC':'F,G5',
         'K_FH':'F,H','B_CC':'G5,G5','K_CH':'G5,H','B_HH':'H,H'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def literal_rational(node):
    """Read only literal arrays/integers/Q fractions; never eval archived code."""
    if isinstance(node,ast.Constant) and isinstance(node.value,(int,str)):
        return node.value
    if isinstance(node,(ast.List,ast.Tuple)):
        return [literal_rational(n) for n in node.elts]
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
        return -literal_rational(node.operand)
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='Q' and not node.keywords:
        return Fraction(*[literal_rational(n) for n in node.args])
    raise ValueError('Unexpected expression in literal dependency binding: '+ast.dump(node))

def supplemental_checks(work):
    fourier=json.loads((work/'second_correction_fourier_check_result.json').read_text())
    endpoint=json.loads((work/'enlarged_endpoint_gate_result.json').read_text())
    constants={}
    for stmt in ast.parse((work/'second_corrected_gate_static.py').read_text()).body:
        if isinstance(stmt,ast.Assign) and len(stmt.targets)==1 and isinstance(stmt.targets[0],ast.Name) and stmt.targets[0].id in ('field_names','caps','fourier_bounds'):
            constants[stmt.targets[0].id]=literal_rational(stmt.value)
    stage12=dict(zip(constants['field_names'],zip(constants['caps'],constants['fourier_bounds'])))
    stage16=dict(zip(endpoint['field_names'],map(Fraction,endpoint['field_critical_caps'])))
    bindings=[]
    for name,pair in CAP_MAP.items():
        bound=Fraction(fourier['pair_fields'][pair]['critical_norm']['rational_upper_bound'])
        cap=Fraction(fourier['coarse_rational_pair_critical_upper_bounds'][pair])
        c12,b12=stage12[name]
        assert bound==b12 and c12==cap==stage16[name] and bound<cap
        bindings.append({'field':name,'pair':pair,'fresh_Fourier_bound':str(bound),'cap':str(cap),
                         'Stage12_literal_bound_matches':True,'Stage12_Stage16_caps_match':True,'strict_cap_pass':True})
    # Full formal bilinear scalar identity, not just rational sample evaluation.
    # Scalar variables in order alpha,beta,gamma,eta,theta. Spatial basis W,F,C,H.
    # The pair symbol X,Y denotes B(X,Y)+B(Y,X) for distinct X,Y and B(X,X)
    # on the diagonal. Exact bilinearity/symmetry alone then determines factors.
    bases={'W':0,'F':1,'G5':2,'H':3}
    def pair(a,b):return ','.join(sorted((a,b),key=bases.get))
    def monomial(i,j):
        e=[0]*5;e[i]+=1;e[j]+=1;return tuple(e)
    direct=defaultdict(int)
    w1=[('W',0),('F',1)]; corr=[('W',2),('G5',3),('H',4)]
    for X,i in w1:
        for Y,j in corr:
            direct[pair(X,Y),monomial(i,j)]+=2 if X==Y else 1
    for index,(X,i) in enumerate(corr):
        for Y,j in corr[index:]:
            direct[pair(X,Y),monomial(i,j)]+=1
    displayed=defaultdict(int)
    for X,Y,i,j,factor in [('W','W',0,2,2),('W','W',2,2,1),('W','F',1,2,1),
                           ('W','G5',0,3,1),('W','G5',2,3,1),('F','G5',1,3,1),
                           ('W','H',0,4,1),('W','H',2,4,1),('F','H',1,4,1),
                           ('G5','G5',3,3,1),('H','H',4,4,1),('G5','H',3,4,1)]:
        displayed[pair(X,Y),monomial(i,j)]+=factor
    assert dict(direct)==dict(displayed) and len(direct)==12
    return {'independent_cross_certificate_cap_binding':bindings,
            'formal_scalar_residual_expansion':'PASS', 'formal_scalar_monomials':len(direct),
            'identity':'C(alpha*W+beta*F,gamma*W+eta*G5+theta*H)+B(correction,correction) equals all nine displayed pair terms',
            'qualification':'Supplement checks literal cap dependencies and bilinear formal scalar expansion; it is not an independent implementation of all Fourier/time-bound algorithms.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    default_source=Path(__file__).resolve().parent.parent/'proof_sources'/'W_research'
    parser.add_argument('--source',type=Path,default=default_source)
    parser.add_argument('--out',type=Path,default=Path('replay_new_W'))
    parser.add_argument('--timeout-per-checker',type=float,default=300)
    args=parser.parse_args();source=args.source.resolve();out=args.out.resolve()
    if source==out or source in out.parents:parser.error('Output must be outside source')
    if out.exists() and any(out.iterdir()):parser.error('Output must be absent or empty')
    out.mkdir(parents=True,exist_ok=True)
    work=out/'working_copy';logs=out/'logs';work.mkdir();logs.mkdir()
    sources=[]
    for directory,script,result in CHECKERS:
        for name in (script,result):
            p=source/directory/name
            assert p.is_file(),directory+'/'+name
            sources.append({'relative_path':str(p.relative_to(source)),'sha256_before':sha(p),'bytes':p.stat().st_size})
        shutil.copyfile(source/directory/script,work/script)
    doc=source/'16_enlarged_endpoint/second_corrected_W_gate.md'
    sources.append({'relative_path':str(doc.relative_to(source)),'sha256_before':sha(doc),'bytes':doc.stat().st_size})
    shutil.copyfile(doc,work/doc.name)
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
    env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0',PYTHONINTMAXSTRDIGITS='0')
    receipt={'scope':'Replay of inspected archived finite exact-arithmetic constructors and archived data-only auditors, in a new isolated copy. Not PDE simulation or a proof assistant.',
             'started_utc':datetime.now(timezone.utc).isoformat(),'source':'source_root','out':'.',
             'path_convention':{'source_root':'Logical root selected by --source; source_files paths are relative to this root.',
                                'out':'Directory containing this receipt; cwd, stdout and stderr paths are relative to it.',
                                'command':'Replayable command in the stated cwd; uses the recorded Python version, not an absolute interpreter path.'},
             'runner_variant':'Portable public wrapper; scientific checker and certificate bytes unchanged. See WRAPPER_CHANGES.txt.',
             'python':sys.version,'platform':platform.platform(),'runner_sha256':sha(Path(__file__)),
             'runtime_environment':{k:env[k] for k in ['PYTHONDONTWRITEBYTECODE','PYTHONHASHSEED','PYTHONINTMAXSTRDIGITS']},
             'original_source_patches':[],'isolated_copy_patches':[],
             'layout_note':'Selected source files flattened byte-for-byte in working_copy; all result JSON files generated afresh, not prepopulated. Paths are basenames in original programs.',
             'expected_unique_checker_count':15,'historical_receipt_count':'10 Picard + 6 endpoint outputs, with one overlapping Fourier output = 15 unique',
             'runs':[],'source_files':sources,
             'limitations':['Replaying the same Fourier constructors does not supply a separate Fourier implementation.',
                            'Archived adversarial programs check encoded data/margins and do not by themselves prove all PDE/functional-analytic lemmas.',
                            'Hardcoded rational field caps are explicitly cross-bound to the fresh Fourier certificate by supplemental_checks.',
                            'Failure of other tested scalar majorants is not nonexistence or divergence.',
                            'Stage16 existence/regularity neighborhoods do not automatically inherit Stage15 heat-seeded Picard convergence.',
                            'No arbitrary smooth 3D datum result, priority claim, external human review or formal proof is established here.']}
    write_json(out/'replay_receipt.json',receipt)
    for number,(directory,script,result) in enumerate(CHECKERS,1):
        command=[sys.executable,str(work/script)]
        target=work/result
        assert not target.exists(),'Output was unexpectedly prepopulated: '+result
        begin=time.monotonic();status={}
        with (logs/(script+'.stdout.txt')).open('wb') as stdout,(logs/(script+'.stderr.txt')).open('wb') as stderr:
            try:
                proc=subprocess.run(command,cwd=work,env=env,stdout=stdout,stderr=stderr,timeout=args.timeout_per_checker)
                exitcode=proc.returncode
            except subprocess.TimeoutExpired:
                exitcode=None;status['timeout']=True
        row={'number':number,'source_script':directory+'/'+script,'command':['python3',script],'cwd':'working_copy',
             'exit_code':exitcode,'runtime_seconds':round(time.monotonic()-begin,6),
             'result':result,'stdout':'logs/'+script+'.stdout.txt','stderr':'logs/'+script+'.stderr.txt',**status}
        if target.exists():
            parsed=json.loads(target.read_text());row.update(result_bytes=target.stat().st_size,
                result_sha256=sha(target),archive_sha256=sha(source/directory/result),
                byte_identical_to_archive=target.read_bytes()==(source/directory/result).read_bytes(),
                result_status=parsed.get('status'))
            for key in ['number_of_assertions','assertion_count','assertions_passed','defect_support_count','mixed_mode_count','mode_count']:
                if key in parsed:row[key]=parsed[key]
        receipt['runs'].append(row);write_json(out/'replay_receipt.json',receipt)
        print(json.dumps({'checker':number,'script':script,'exit_code':exitcode,'runtime_seconds':row['runtime_seconds'],'status':row.get('result_status'),'byte_identical':row.get('byte_identical_to_archive')}),flush=True)
        if exitcode!=0 or not row.get('byte_identical_to_archive'):
            receipt['status']='STOPPED_CHECKER_FAILURE_OR_OUTPUT_MISMATCH';write_json(out/'replay_receipt.json',receipt)
            return 1
    receipt['supplemental_checks']=supplemental_checks(work)
    for row in receipt['source_files']:
        row['sha256_after']=sha(source/row['relative_path'])
        row['unchanged']=row['sha256_before']==row['sha256_after']
    receipt.update(status='PASS_ALL_15_REPLAY_OUTPUTS_BYTE_IDENTICAL',
                   checker_count=len(receipt['runs']),all_selected_originals_unchanged=all(r['unchanged'] for r in receipt['source_files']),
                   all_result_bytes_identical=all(r['byte_identical_to_archive'] for r in receipt['runs']),
                   total_checker_runtime_seconds=sum(r['runtime_seconds'] for r in receipt['runs']),
                   completed_utc=datetime.now(timezone.utc).isoformat())
    write_json(out/'replay_receipt.json',receipt)
    print(json.dumps({'status':receipt['status'],'runtime_seconds':receipt['total_checker_runtime_seconds'],'receipt':'replay_receipt.json'}),flush=True)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
