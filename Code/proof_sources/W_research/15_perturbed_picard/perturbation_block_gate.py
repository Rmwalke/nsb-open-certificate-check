"""Fresh exact block perturbation certificate; historical JSON data only."""
from fractions import Fraction as F
from math import factorial, isqrt
from functools import lru_cache
from pathlib import Path
import hashlib
import json

SOURCE=Path("refined_W_block_gate_result.json")
RAW=SOURCE.read_bytes()
EXPECTED="17b578bd37cb971cfbd21a58749bc95a94d62e126fcced0fe21bb3b85a3420b0"
assert hashlib.sha256(RAW).hexdigest()==EXPECTED
DATA=json.loads(RAW)
PASS=[t for t in DATA["tests"] if t["passed"]]
assert len(PASS)==1 and PASS[0]["reference"]=="w3"
BASE=PASS[0]
ROWS=BASE["rows"]
assert len(ROWS)==17
C=F(707106781187,2000000000000)
assert C*C>F(1,8)
S=10**18
ES=10**60

def down(x,scale=S):
    return F(x.numerator*scale//x.denominator,scale)
def up(x,scale=S):
    return F((x.numerator*scale+x.denominator-1)//x.denominator,scale)
def root(x):
    assert x>=0
    n=isqrt(x.numerator*S*S//x.denominator)
    lo=F(n,S)
    return lo,lo if lo*lo==x else F(n+1,S)

@lru_cache(None)
def exp_interval(x):
    assert x>=0
    n=0
    while x>1:
        x/=2
        n+=1
    hi=sum(((-x)**j/F(factorial(j)) for j in range(65)),F(0))
    lo=hi+(-x)**65/F(factorial(65))
    lo=down(lo,ES); hi=up(hi,ES)
    assert 0<lo<=hi
    for _ in range(n):
        lo=down(lo*lo,ES); hi=up(hi*hi,ES)
    return lo,hi

@lru_cache(None)
def heat_factor(a,length):
    if a==0:
        return F(1)
    best=exp_interval(a)[1]
    if length is None:
        return best
    q=length/(2*a+length)
    loglo=2*sum((q**(2*j+1)/F(2*j+1) for j in range(24)),F(0))
    loghi=loglo+2*q**49/(49*(1-q*q))
    all_lambda=length/(a+length)*exp_interval((a/length)*loglo)[1]
    best=min(best,root(all_lambda)[1])
    if loghi<=2*length:
        boundary=exp_interval(2*a)[1]*(1-exp_interval(2*length)[0])
        best=min(best,root(boundary)[1])
    return best

VS=[F(r["V_upper"]) for r in ROWS]
DS=[F(r["defect_upper"]) for r in ROWS]
KS=[[F(x) for x in r["coefficients"]]+[F(1)] for r in ROWS]
GRID=[F(x) if x!="infinity" else None for x in BASE["grid"]]
INITIAL_CAPS=[2*C*heat_factor(a,None if b is None else b-a)
              for a,b in zip(GRID[:-1],GRID[1:])]
assert all(0<x<=2*C for x in INITIAL_CAPS)

# Recheck the imported baseline invariant-set algebra before using its data.
for i,r in enumerate(BASE["radii_upper"]):
    radius=F(r)
    old=[F(x) for x in BASE["radii_upper"][:i+1]]
    assert DS[i]+C*sum((KS[i][j]*(2*VS[j]*old[j]+old[j]**2)
                       for j in range(i+1)),F(0))<=radius
    assert 2*C*(VS[i]+radius)<F(3,5)

def certify(chi,keep_rows=True):
    low=[]; high=[]; report=[]
    for i in range(17):
        d=DS[i]+chi*INITIAL_CAPS[i]
        hl=down(C*sum((KS[i][j]*(2*VS[j]*low[j]+low[j]**2)
                      for j in range(i)),F(0)))
        hh=up(C*sum((KS[i][j]*(2*VS[j]*high[j]+high[j]**2)
                    for j in range(i)),F(0)))
        gamma=1-2*C*VS[i]
        disclo=gamma*gamma-4*C*(d+hh)
        dischi=gamma*gamma-4*C*(d+hl)
        row={"block":i+1,"perturbed_defect_cap":str(d),
             "initial_heat_cap":str(chi*INITIAL_CAPS[i]),
             "history_lower":str(hl),"history_upper":str(hh),
             "discriminant_lower":str(disclo),"discriminant_upper":str(dischi)}
        if gamma<=0 or dischi<=0:
            row["status"]="DECLARED_MAJORANT_FAILS"
            report.append(row)
            return {"chi":str(chi),"passed":False,"failure_proved":True,
                    "first_failed_block":i+1,"rows":report if keep_rows else [],
                    "failed_discriminant_upper":str(dischi)}
        assert disclo>0
        rl=(gamma-root(dischi)[1])/(2*C)
        rh=(gamma-root(disclo)[0])/(2*C)
        assert 0<=rl<=rh
        assert d+hh+2*C*VS[i]*rh+C*rh*rh<=rh
        assert chi*INITIAL_CAPS[i]<=rh
        assert 2*C*(VS[i]+rh)<1
        low.append(rl); high.append(rh)
        row.update({"status":"PASS","radius_lower":str(rl),"radius_upper":str(rh)})
        report.append(row)
    matrix=[[2*C*KS[i][j]*(VS[j]+high[j]) for j in range(i+1)]
            for i in range(17)]
    diag=max(matrix[i][i] for i in range(17))
    budget=sum(((VS[i]+high[i])**2 for i in range(17)),F(0))
    return {"chi":str(chi),"passed":True,"rows":report if keep_rows else [],
            "radii_upper":list(map(str,high)),
            "contraction_matrix":[list(map(str,row)) for row in matrix],
            "spectral_radius_upper":str(diag),
            "solution_X_squared_upper":str(budget)}

TRIALS=[F(0),F(1,100000),F(1,50000),F(1,10000),F(1,5000),
        F(1,2000),F(1,1000),F(1,500),F(1,200),F(1,100),F(1,20)]
tests=[]
for chi in TRIALS:
    test=certify(chi)
    tests.append(test)
    print(json.dumps({"chi":str(chi),"passed":test["passed"],
        "diagonal":float(F(test["spectral_radius_upper"])) if test["passed"] else None,
        "budget":float(F(test["solution_X_squared_upper"])) if test["passed"] else None,
        "failed_block":test.get("first_failed_block")}),flush=True)
# Use a robust predeclared margin, rather than the largest merely passing radius.
eligible=[t for t in tests if t["passed"] and F(t["chi"])>0
          and F(t["spectral_radius_upper"])<F(2,3)
          and F(t["solution_X_squared_upper"])<12]
assert eligible
selected=max(eligible,key=lambda t:F(t["chi"]))
matrix=[[F(x) for x in row] for row in selected["contraction_matrix"]]
radii=[F(x) for x in selected["radii_upper"]]
q=F(3,4)
p=[]
for i,row in enumerate(matrix):
    assert row[i]<q
    p.append((1+sum((row[j]*p[j] for j in range(i)),F(0)))/(q-row[i]))
    assert sum((row[j]*p[j] for j in range(i+1)),F(0))==q*p[i]-1
alpha=max(r/p_i for r,p_i in zip(radii,p))
alpha_shift=max((r+F(selected["chi"])*v)/p_i
                for r,v,p_i in zip(radii,INITIAL_CAPS,p))
constant2=alpha*alpha*sum((x*x for x in p),F(0))
constant_shift2=alpha_shift*alpha_shift*sum((x*x for x in p),F(0))
integer_constant=root(constant2)[1].__ceil__()
integer_shift=root(constant_shift2)[1].__ceil__()
assert constant2<=integer_constant**2
assert constant_shift2<=integer_shift**2
rate={"q":str(q),"inverse_weights":list(map(str,p)),
      "scale_from_w3":str(alpha),"constant_squared_from_w3":str(constant2),
      "integer_constant_from_w3":integer_constant,
      "scale_from_w3_plus_Sh":str(alpha_shift),
      "constant_squared_from_w3_plus_Sh":str(constant_shift2),
      "integer_constant_from_w3_plus_Sh":integer_shift}
result={"status":"PASS_PERTURBED_W3_BLOCK_CONVERGENCE",
        "source_json":str(SOURCE),"source_sha256":EXPECTED,
        "scope":"Real, mean-zero, divergence-free smooth h; chi(h)=A_-1(h).",
        "c_upper":str(C),"grid":BASE["grid"],
        "reference_block_caps":list(map(str,VS)),
        "baseline_defect_block_caps":list(map(str,DS)),
        "source_memory_factors":[list(map(str,row)) for row in KS],
        "heat_profile_caps":[str(x) for x in INITIAL_CAPS],
        "selected_chi":selected["chi"],"selected_certificate":selected,
        "rate":rate,"bounded_trials":tests,
        "not_claimed":["No enlargement over the older chi<=1/20 regularity tube.",
          "Standard heat-seeded Picard entry requires a separate finite-prefix certificate.",
          "Failure of a declared majorant does not imply divergence or failure of other gates."]}
Path("perturbation_block_gate_result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"chi":selected["chi"],
                  "rate":str(q),"constant_from_w3":integer_constant,
                  "constant_from_shifted_reference":integer_shift}),flush=True)

