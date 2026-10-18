"""Fresh exact endpoint-moment checks for w2; no historical source imports."""
from fractions import Fraction as F
from math import factorial, isqrt
from functools import lru_cache
from pathlib import Path
import hashlib
import json

SOURCE = Path("third_heat_defect_static_result.json")
INPUT = json.loads(SOURCE.read_text())
S = 10**18
ES = 10**60
CBAR = F(707106781187,2000000000000)
assert CBAR*CBAR > F(1,8)

def down(q,scale=S):
    return F((q.numerator*scale)//q.denominator,scale)
def up(q,scale=S):
    return F((q.numerator*scale+q.denominator-1)//q.denominator,scale)
def sqrt_interval(q):
    assert q >= 0
    n = isqrt((q.numerator*S*S)//q.denominator)
    lo=F(n,S)
    return lo, lo if lo*lo==q else F(n+1,S)

@lru_cache(None)
def exp_minus(q):
    assert q >= 0
    count=0
    while q>1:
        q/=2
        count+=1
    even=sum(((-q)**n/F(factorial(n)) for n in range(65)),F(0))
    odd=even+(-q)**65/F(factorial(65))
    lo=down(odd,ES); hi=up(even,ES)
    assert 0<lo<=hi
    for _ in range(count):
        lo=down(lo*lo,ES); hi=up(hi*hi,ES)
    return lo,hi

def add(p,q):
    r=dict(p)
    for key,val in q.items():
        r[key]=r.get(key,F(0))+val
        if not r[key]: del r[key]
    return r

def scaled(p,c):
    return {k:c*v for k,v in p.items() if c*v}

def product(p,q):
    r={}
    for (a,d),x in p.items():
        for (b,e),y in q.items():
            key=(a+b,d+e)
            r[key]=r.get(key,F(0))+x*y
    return {k:v for k,v in r.items() if v}

def total_integral(p):
    return sum((v*F(factorial(d),a**(d+1))
                for (a,d),v in p.items()),F(0))

def tail_poly(p):
    result={}
    for (rate,degree),value in p.items():
        for power in range(degree+1):
            key=(rate,power)
            result[key]=result.get(key,F(0))+value*F(
                factorial(degree),factorial(power)*rate**(degree-power+1))
    return {k:v for k,v in result.items() if v}

def frozen(p):
    return tuple(sorted((a,d,v) for (a,d),v in p.items()))

@lru_cache(None)
def tail_interval(frozen_tail,t):
    if t is None: return F(0),F(0)
    by_rate={}
    for a,d,v in frozen_tail:
        by_rate[a]=by_rate.get(a,F(0))+v*t**d
    lo=hi=F(0)
    for a,v in by_rate.items():
        el,eh=exp_minus(a*t)
        if v>=0: lo+=v*el; hi+=v*eh
        else: lo+=v*eh; hi+=v*el
    return lo,hi

def block_integral(frozen_tail,left,right):
    x=tail_interval(frozen_tail,left)
    y=tail_interval(frozen_tail,right)
    return x[0]-y[1],x[1]-y[0]

# Reconstruct every signed Fourier-mode squared norm directly from JSON data.
groups={}
for entry in INPUT["all_defect_modes"]:
    square={}
    for component in entry["imaginary_vector_profiles"]:
        p={(int(t["rate"]),int(t["power"])):F(t["coefficient"])
           for t in component}
        square=add(square,product(p,p))
    assert total_integral(square)==F(entry["exact_L2_time_squared"])
    key=frozen(tail_poly(square))
    groups[key]=groups.get(key,0)+1
assert sum(groups.values())==342

SQRT3LO,SQRT3HI=sqrt_interval(F(3))
assert SQRT3LO**2<=3<=SQRT3HI**2
L6={(4,2):F(1,2),(4,1):F(-1,2),(4,0):F(1,4),(6,0):F(-1,4)}
base={(1,0):F(21,8),(2,1):F(3),(3,1):F(3,2),(5,0):F(3,8)}
Ulo=add(base,scaled(L6,SQRT3LO))
Uhi=add(base,scaled(L6,SQRT3HI))
Ulo2=product(Ulo,Ulo); Uhi2=product(Uhi,Uhi)
UTLO=frozen(tail_poly(Ulo2)); UTHI=frozen(tail_poly(Uhi2))
PLO=total_integral(Ulo2); PHI=total_integral(Uhi2)
assert F(441,128)<PLO<=PHI<F(217,32)

@lru_cache(None)
def block_data(a,b):
    vlo=sqrt_interval(max(F(0),block_integral(UTLO,a,b)[0]))[0]
    vhi=sqrt_interval(block_integral(UTHI,a,b)[1])[1]
    dhi=F(0)
    for p,count in groups.items():
        _,upper=block_integral(p,a,b)
        assert upper>=0
        dhi+=count*sqrt_interval(upper)[1]
    return vlo,vhi,dhi

@lru_cache(None)
def kernel(gap,duration):
    if gap==0: return F(1)
    result=exp_minus(gap)[1]
    if duration is None: return result
    q=duration/(2*gap+duration)
    nterms=24
    loglo=2*sum((q**(2*n+1)/F(2*n+1)
                 for n in range(nterms)),F(0))
    loghi=loglo+2*q**(2*nterms+1)/(
        F(2*nterms+1)*(1-q*q))
    unconstrained=(duration/(gap+duration)
                   *exp_minus((gap/duration)*loglo)[1])
    result=min(result,sqrt_interval(unconstrained)[1])
    if loghi<=2*duration:
        constrained=(exp_minus(2*gap)[1]
                     *(1-exp_minus(2*duration)[0]))
        result=min(result,sqrt_interval(constrained)[1])
    return result

GRIDS=[
 ("prior17_exact_U2",[F(0),F(1,16),F(1,8),F(3,16),F(1,4),
  F(5,16),F(3,8),F(7,16),F(1,2),F(5,8),F(3,4),F(7,8),
  F(1),F(5,4),F(3,2),F(2),F(3),None]),
 ("fine43_exact_U2",[F(j,16) for j in range(33)]
  +[F(j,8) for j in range(17,25)]+[F(4),F(6),None]),
 ("medium17_exact_U2",[F(j,8) for j in range(9)]
  +[F(j,4) for j in range(5,9)]+[F(5,2),F(3),F(4),F(6),None]),
]

def test_grid(name,grid,refined_diagonal=False,fourth_data=None):
    rows=[]; lows=[]; highs=[]; Vs=[]; Ks=[]; diagonals=[]
    for i,(a,b) in enumerate(zip(grid[:-1],grid[1:])):
        vl,V,d=block_data(a,b)
        reference2_upper=V
        third_defect_upper=d
        if fourth_data is not None:
            source_row=fourth_data["block_defect_bounds"][i]
            assert F(source_row["left"])==a
            if b is not None: assert F(source_row["right"])==b
            vl=max(F(0),vl-d)
            V=V+d
            d=F(source_row["full_defect_norm_upper"])
        Vs.append(V)
        coeffs=[kernel(a-grid[j+1],None if b is None else b-a)
                for j in range(i)]
        hlo=down(CBAR*sum((k*(2*Vs[j]*lows[j]+lows[j]**2)
                          for j,k in enumerate(coeffs)),F(0)))
        hhi=up(CBAR*sum((k*(2*Vs[j]*highs[j]+highs[j]**2)
                        for j,k in enumerate(coeffs)),F(0)))
        diagonal=2*CBAR*V
        M_upper=None
        if refined_diagonal and b is not None:
            M_upper=3*exp_minus(a)[1]+3*a*exp_minus(2*a)[1]
            diagonal=min(diagonal,F(17,25)*M_upper*sqrt_interval(b-a)[1])
        diagonals.append(diagonal)
        gamma=1-diagonal
        disc_lo=gamma*gamma-4*CBAR*(d+hhi)
        disc_hi=gamma*gamma-4*CBAR*(d+hlo)
        row={"block":i+1,"left":str(a),
             "right":"infinity" if b is None else str(b),
             "reference2_norm_upper":str(reference2_upper),
             "third_heat_defect_upper":str(third_defect_upper),
             "V_lower":str(vl),"V_upper":str(V),
             "defect_upper":str(d),"history_lower":str(hlo),
             "history_upper":str(hhi),"gamma":str(gamma),
             "diagonal_linear_upper":str(diagonal),
             "sup_U1_upper":None if M_upper is None else str(M_upper),
             "discriminant_lower":str(disc_lo),
             "discriminant_upper":str(disc_hi),
             "coefficients":list(map(str,coeffs))}
        if gamma<=0 or disc_hi<=0:
            row["status"]="DECLARED_MAJORANT_PROVABLY_FAILS"
            rows.append(row)
            return {"name":name,"reference":"w3" if fourth_data else "w2",
                    "refined_diagonal":refined_diagonal,
                    "passed":False,"failure_proved":True,
                    "first_failed_block":i+1,"rows":rows,
                    "grid":[str(x) if x is not None else "infinity" for x in grid]}
        assert disc_lo>0, "Insufficient interval precision; no false pass allowed"
        rl=(gamma-sqrt_interval(disc_hi)[1])/(2*CBAR)
        rh=(gamma-sqrt_interval(disc_lo)[0])/(2*CBAR)
        assert 0<=rl<=rh
        assert d+hhi+diagonal*rh+CBAR*rh*rh<=rh
        assert diagonal+2*CBAR*rh<1
        lows.append(rl); highs.append(rh); Ks.append(coeffs+[F(1)])
        row.update({"status":"PASS","radius_lower":str(rl),"radius_upper":str(rh)})
        rows.append(row)
    matrix=[[(diagonals[i]+2*CBAR*highs[i] if i==j else
              2*CBAR*Ks[i][j]*(Vs[j]+highs[j]))
             for j in range(i+1)] for i in range(len(Vs))]
    assert max(matrix[i][i] for i in range(len(Vs)))<1
    solution_budget=sum(((V+r)**2 for V,r in zip(Vs,highs)),F(0))
    if fourth_data is not None:
        assert max(matrix[i][i] for i in range(len(Vs)))<F(3,5)
        assert solution_budget<8
    return {"name":name,"reference":"w3" if fourth_data else "w2",
            "refined_diagonal":refined_diagonal,
            "passed":True,"rows":rows,
            "grid":[str(x) if x is not None else "infinity" for x in grid],
            "radii_upper":list(map(str,highs)),
            "contraction_matrix":[list(map(str,row)) for row in matrix],
            "solution_X_squared_upper":str(solution_budget),
            "spectral_radius_upper":str(max(matrix[i][i] for i in range(len(Vs))))}

result={"scope":"Fresh rational endpoint moments; no evolution or quadrature.",
        "source_json":str(SOURCE),
        "source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "mode_count":342,"distinct_squared_norm_profiles":len(groups),
        "sqrt3_interval":[str(SQRT3LO),str(SQRT3HI)],
        "P2_exact_enclosure":[str(PLO),str(PHI)],"tests":[]}
for refinement in (False,True):
    for name,grid in GRIDS:
        test=test_grid(name,grid,refinement)
        result["tests"].append(test)
        last=test["rows"][-1]
        print(json.dumps({"grid":name,"refined_diagonal":refinement,
                          "passed":test["passed"],
                          "first_failed_block":test.get("first_failed_block"),
                          "last_interval":[last["left"],last["right"]],
                          "last_discriminant_upper_for_orientation":
                           float(F(last["discriminant_upper"]))}),flush=True)
result["convergence_certified"]=any(t["passed"] for t in result["tests"])
fourth_source=Path("fourth_heat_defect_static_result.json")
if fourth_source.exists():
    fourth_bytes=fourth_source.read_bytes()
    fourth_data=json.loads(fourth_bytes)
    result["fourth_source_sha256"]=hashlib.sha256(fourth_bytes).hexdigest()
    name,grid=GRIDS[0]
    test=test_grid("w3_prior17_exact_U2_plus_D3",grid,False,fourth_data)
    result["tests"].append(test)
    last=test["rows"][-1]
    print(json.dumps({"grid":test["name"],"passed":test["passed"],
                      "first_failed_block":test.get("first_failed_block"),
                      "last_interval":[last["left"],last["right"]],
                      "last_discriminant_upper_for_orientation":
                       float(F(last["discriminant_upper"]))}),flush=True)
result["convergence_certified"]=any(t["passed"] for t in result["tests"])
result["status"]=("PASS_REFINED_CONVERGENCE_CERTIFICATE" if
                  result["convergence_certified"] else
                  "PASS_REFINED_FAILURE_DIAGNOSTICS")
Path("refined_W_block_gate_result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],
                  "P2_interval_for_orientation":[float(PLO),float(PHI)],
                  "groups":len(groups)}),flush=True)
