"""Independent new exact-algebra audit of enlarged endpoint gate data.

Reads frozen JSON and source hashes only. No historical research code is
executed or imported; no time quadrature or evolution is computed.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json


def plus(*args):
    p = {}
    for q in args:
        for key,v in q.items():
            p[key] = p.get(key,F(0))+v
    return {key:v for key,v in p.items() if v}


def times_scalar(c,p):
    return {key:c*v for key,v in p.items() if c*v}


def product(p,q):
    result = {}
    for (rate,degree),v in p.items():
        for (rate2,degree2),w in q.items():
            key = (rate+rate2,degree+degree2)
            result[key] = result.get(key,F(0))+v*w
    return {key:v for key,v in result.items() if v}


def integral(p,n=0):
    return sum((v*F(factorial(d),(r+2*n)**(d+1)) for (r,d),v in p.items()),F(0))


def exponential_upper(x):
    assert 0 <= x < 102
    terms = [x**j/F(factorial(j)) for j in range(101)]
    return sum(terms,F(0)) + terms[-1]*x/101/(1-x/102)


def primitive_upper_direct(s):
    """Independent direct atanh series; no dyadic normalization."""
    assert 0<=s<2
    q=s/(4-s)
    n=300
    series=sum((q**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    remainder=q**(2*n+1)/((2*n+1)*(1-q*q))
    return -2*s+8*(series+remainder)


def main():
    name = 'enlarged_endpoint_gate_result.json'
    raw = Path(name).read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    assert sha == '661a9c1c426b59e5f9d0a9477628f03f3b96a1ec8f86dc69f4c7d5e919ae4f0f'
    data = json.loads(raw)
    for source, expected in data['source_sha256'].items():
        assert hashlib.sha256(Path(source).read_bytes()).hexdigest() == expected
    caps = list(map(F,data['field_critical_caps']))
    assert caps == list(map(F,['17/8','171/100','27/20','87/40','6/5','183/100','11/25','3/4','9/25']))
    lower,upper = map(F,data['sqrt3_interval'])
    assert 0 < lower and lower*lower < 3 < upper*upper
    l1 = {(1,0):F(1,4),(3,0):F(-1,4),(3,1):F(-1,2)}
    l5 = {(3,0):F(-1,4),(3,1):F(1,2),(5,0):F(1,4)}
    l6 = {(4,0):F(1,4),(4,1):F(-1,2),(4,2):F(1,2),(6,0):F(-1,4)}
    verified = {}
    for record in data['profiles']:
        m = F(record['amplitude'])
        assert 0 < m and m*m < F(8,3) and m*m < 3
        x={(1,0):m}; y={(2,1):m*m}
        p=times_scalar(m**3/2,l1); q=times_scalar(m**3,l5); r=times_scalar(m**4,l6)
        diff=plus(x,times_scalar(-1,p))
        terms=[plus(times_scalar(2,product(x,p)),times_scalar(-1,product(p,p))),
               product(diff,q),product(diff,r),product(y,p),product(y,q),
               product(y,r),product(q,q),product(q,r),product(r,r)]
        residual=plus(*(times_scalar(c,v) for c,v in zip(caps,terms)))
        saved={(v['rate'],v['degree']):F(v['coefficient']) for v in record['residual_profile']}
        assert residual == saved
        reference=plus(times_scalar(3,diff),times_scalar(3,y),times_scalar(F(3,2),q))
        lo=plus(reference,times_scalar(lower,r)); hi=plus(reference,times_scalar(upper,r))
        assert integral(product(lo,lo)) == F(record['P2_lower'])
        P=integral(product(hi,hi))
        assert P == F(record['P2_upper']) < F(record['P1_upper'])
        assert integral(residual) == F(record['residual_total_upper'])
        c=F(9,8)*m*m; delta=m**3/2+F(9,128)*m**4
        assert c==F(record['c']) and delta==F(record['delta'])
        moments=[integral(residual,n) for n in range(26)]
        assert all(v>0 for v in moments)
        assert all(moments[j+1]<=moments[j] for j in range(25))
        series=sum((c**n/F(factorial(n))*moments[n] for n in range(25)),F(0))
        series+=c**25/F(factorial(25))*moments[25]/(1-c/26)
        assert series==F(record['positive_series_upper'])
        assert exponential_upper(delta)*series<=F(record['weighted_residual_upper'])
        verified[m]=record
    passing=[]
    for record in data['finite_candidate_classes']:
        m,radius=F(record['amplitude']),F(record['chi'])
        profile=verified[m]
        P=F(profile['P2_upper'])
        assert 0<=radius<2
        assert F(record['P2_upper'])==P
        assert F(record['residual_total_upper'])==F(profile['residual_total_upper'])
        assert F(record['weighted_residual_upper'])==F(profile['weighted_residual_upper'])
        initial=radius*exponential_upper(P/4-radius/2)
        assert initial<=F(record['initial_endpoint_upper'])
        M=F(record['endpoint_upper'])
        assert M==F(record['initial_endpoint_upper'])+F(record['weighted_residual_upper'])
        if record['passed']:
            B=F(record['selected_barrier'])
            assert radius<B<2
            safe_psi=B/exponential_upper(B/2)
            assert F(record['Psi_barrier_lower'])<=safe_psi
            assert M<F(record['Psi_barrier_lower'])
            assert F(record['endpoint_margin_lower'])==F(record['Psi_barrier_lower'])-M
            passing.append({'amplitude':str(m),'chi':str(radius),'barrier':str(B),
                            'endpoint_upper_orientation_only':float(M),
                            'strict_gate_margin_orientation_only':float(F(record['endpoint_margin_lower']))})
    assert any(v['amplitude']=='1' and v['chi']=='1/10' for v in passing)
    assert any(v['amplitude']=='8/7' and v['chi']=='1/100' for v in passing)
    budget_name='enlarged_budget_static_result.json'
    budget_raw=Path(budget_name).read_bytes()
    budget_data=json.loads(budget_raw)
    assert budget_data['source_sha256']==sha
    budget_rows=budget_data['new_classes']
    assert len(budget_rows)==5
    original_by_pair={(F(v['amplitude']),F(v['chi'])):v for v in data['finite_candidate_classes'] if v['passed']}
    verified_budgets=[]
    for row in budget_rows+[budget_data['old_classes_common_corollary']]:
        radius=F(row['chi_radius']); B=F(row['critical_barrier'])
        P=F(row['reference_budget_upper']); R=F(row['residual_integral_upper'])
        delta=1-B/2
        assert 0<=radius<B<2 and delta==F(row['delta'])
        fb=F(row['F_barrier_upper']); fr=F(row['F_initial_upper'])
        assert fb>=primitive_upper_direct(B)>=0
        assert fr>=primitive_upper_direct(radius)>=0
        C=fr+B*R/delta
        assert C==F(row['forcing_primitive_cost_upper'])
        rootQ=F(row['error_X_norm_upper'])
        d=rootQ*rootQ-(2*fb+C)
        assert rootQ>=0 and d>=0 and d*d>=4*fb*(fb+C)
        G=F(row['solution_X_squared_upper'])
        d=G-P-rootQ*rootQ
        assert d>=0 and d*d>=4*P*rootQ*rootQ
        exponent=row['integer_NSB_exponent']
        assert G<exponent and exponent>=2
        if row['label']=='new endpoint class':
            old=original_by_pair[(F(row['amplitude']),radius)]
            assert B==F(old['selected_barrier'])
            assert P==F(old['P2_upper']) and R==F(old['residual_total_upper'])
        else:
            assert B==1 and radius==F(1,20) and P==F(213,25) and R==F(3,20)
            assert exponent==24
        verified_budgets.append({'amplitude':row['amplitude'],'chi_radius':str(radius),
                                'barrier':str(B),'NSB_exponent':exponent,
                                'budget_upper_orientation_only':float(G)})
    result={'status':'PASS_INDEPENDENT_ENLARGED_ENDPOINT_AUDIT',
            'source_json':name,'source_sha256':sha,
            'profiles_independently_rebuilt':len(verified),
            'passing_classes':passing,
            'budget_source_json':budget_name,
            'budget_source_sha256':hashlib.sha256(budget_raw).hexdigest(),
            'independently_verified_budgets':verified_budgets,
            'new_budget_lemma':'Independently proved in enlarged_endpoint_adversarial.md; direct rational log and radical checks PASS.',
            'scope':'Smooth real mean-zero divergence-free aW+h. Separate rows form a union, not a merged rectangle.'}
    Path('enlarged_endpoint_adversarial_static_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
