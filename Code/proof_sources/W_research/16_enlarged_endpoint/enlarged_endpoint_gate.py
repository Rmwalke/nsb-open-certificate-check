"""New exact endpoint-class tests; no historical code imported or executed."""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
from functools import lru_cache
import hashlib
import json

C_NAMES=["F","K_WC","K_WH","G","K_FC","K_FH","B_CC","K_CH","B_HH"]
CAPS=[F(17,8),F(171,100),F(27,20),F(87,40),F(6,5),
      F(183,100),F(11,25),F(3,4),F(9,25)]
L1={(1,0):F(1,4),(3,0):F(-1,4),(3,1):F(-1,2)}
L5={(5,0):F(1,4),(3,0):F(-1,4),(3,1):F(1,2)}
L6={(4,2):F(1,2),(4,1):F(-1,2),(4,0):F(1,4),(6,0):F(-1,4)}
S=10**18
n=isqrt(3*S*S)
SQRT3LO=F(n,S); SQRT3HI=F(n+1,S)
assert SQRT3LO**2<3<SQRT3HI**2
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def add(*polys):
    out={}
    for p in polys:
        for key,x in p.items():out[key]=out.get(key,F(0))+x
    return {k:v for k,v in out.items() if v}
def scale(p,x):return {k:v*x for k,v in p.items() if v*x}
def multiply(p,q):
    out={}
    for (a,m),x in p.items():
        for (b,n),y in q.items():
            k=(a+b,m+n);out[k]=out.get(k,F(0))+x*y
    return {k:v for k,v in out.items() if v}
def moment(p,n=0):
    return sum((x*F(factorial(d),(a+2*n)**(d+1))
                for (a,d),x in p.items()),F(0))
def exp_upper(x):
    check("nonnegative exponential argument "+str(x),x>=0)
    term=total=F(1)
    for n in range(1,65):
        term=term*x/n
        total+=term
    check("exponential geometric tail ratio "+str(x),x<66)
    return total+(term*x/65)/(1-x/66)
def psi_lower(b):return b/exp_upper(b/2)

@lru_cache(None)
def profile_certificate(m):
    check("U2 amplitude monotonicity "+str(m),m*m<F(8,3))
    check("residual amplitude monotonicity "+str(m),m*m<3)
    x={(1,0):m}; y={(2,1):m*m}
    p=scale(L1,m**3/2); q=scale(L5,m**3); r=scale(L6,m**4)
    xm=add(x,scale(p,-1))
    profiles=[
        add(scale(multiply(x,p),2),scale(multiply(p,p),-1)),
        multiply(xm,q),multiply(xm,r),multiply(y,p),multiply(y,q),
        multiply(y,r),multiply(q,q),multiply(q,r),multiply(r,r)]
    R=add(*(scale(p,c) for p,c in zip(profiles,CAPS)))
    Ubase=add(scale(xm,3),scale(y,3),scale(q,F(3,2)))
    Ulo=add(Ubase,scale(r,SQRT3LO))
    Uhi=add(Ubase,scale(r,SQRT3HI))
    Plo=moment(multiply(Ulo,Ulo)); Phi=moment(multiply(Uhi,Uhi))
    P1=F(9,2)*m*m+2*m**3+F(9,32)*m**4
    check("P2 bracket "+str(m),0<Plo<=Phi<P1)
    c=F(9,8)*m*m
    delta=m**3/2+F(9,128)*m**4
    terms=[]; coeff=F(1); prev=None
    for n in range(26):
        if n:coeff=coeff*c/n
        I=moment(R,n)
        check("positive residual moment "+str((m,n)),I>0)
        if prev is not None:check("decreasing residual moment "+str((m,n)),I<=prev)
        prev=I;terms.append(coeff*I)
    check("positive-series tail ratio "+str(m),c<26)
    series=sum(terms[:25],F(0))+terms[25]/(1-c/26)
    Mres=exp_upper(delta)*series
    return {"amplitude":str(m),"P2_lower":str(Plo),"P2_upper":str(Phi),
      "P1_upper":str(P1),"residual_total_upper":str(moment(R)),
      "c":str(c),"delta":str(delta),"positive_series_upper":str(series),
      "weighted_residual_upper":str(Mres),
      "residual_profile":[{"rate":a,"degree":d,"coefficient":str(v)}
           for (a,d),v in sorted(R.items())]}

CANDIDATES=[
    (F(1),F(2,25)),(F(1),F(1,10)),(F(1),F(11,100)),
    (F(11,10),F(1,40)),(F(11,10),F(1,20)),
    (F(9,8),F(1,100)),(F(9,8),F(1,40)),
    (F(8,7),F(1,100)),
    (F(7,6),F(1,100)),(F(6,5),F(1,100)),(F(5,4),F(0))]
BARRIERS=[F(1),F(11,10),F(6,5),F(5,4),F(7,5),F(3,2),
          F(8,5),F(7,4),F(9,5),F(19,10)]
PSI_CAPS={str(B):psi_lower(B) for B in BARRIERS}
records=[]
for m,chi in CANDIDATES:
    base=profile_certificate(m)
    P=F(base["P2_upper"]); MR=F(base["weighted_residual_upper"])
    check("initial perturbation increasing branch "+str((m,chi)),0<=chi<2)
    arg=P/4-chi/2
    check("positive initial exponential exponent "+str((m,chi)),arg>=0)
    initial=chi*exp_upper(arg)
    M=initial+MR
    passes=[B for B in BARRIERS if chi<B and M<PSI_CAPS[str(B)]]
    record={"amplitude":str(m),"chi":str(chi),"passed":bool(passes),
       "P2_upper":str(P),"residual_total_upper":base["residual_total_upper"],
       "weighted_residual_upper":str(MR),"initial_endpoint_upper":str(initial),
       "endpoint_upper":str(M),"passing_barriers":list(map(str,passes))}
    if passes:
        B=passes[0]
        record.update({"selected_barrier":str(B),
                       "Psi_barrier_lower":str(PSI_CAPS[str(B)]),
                       "endpoint_margin_lower":str(PSI_CAPS[str(B)]-M)})
    records.append(record)
    print(json.dumps({"m":str(m),"chi":str(chi),"pass":bool(passes),
         "B":str(passes[0]) if passes else None,
         "M":float(M),"MR":float(MR),"P":float(P)}),flush=True)

for m,chi in [(F(1),F(1,10)),(F(9,8),F(1,100))]:
    found=[x for x in records if F(x["amplitude"])==m and F(x["chi"])==chi]
    check("target enlargement "+str((m,chi)),len(found)==1 and found[0]["passed"])

sources=["second_corrected_W_gate.md","second_corrected_gate_static_result.json",
         "second_correction_fourier_check_result.json"]
result={"status":"PASS_ENLARGED_ENDPOINT_REGULARITY_CLASSES",
        "scope":"Smooth real mean-zero divergence-free u0=aW+h, viscosity1, full Fourier support h.",
        "sqrt3_interval":[str(SQRT3LO),str(SQRT3HI)],
        "field_names":C_NAMES,"field_critical_caps":list(map(str,CAPS)),
        "source_sha256":{name:hashlib.sha256(Path(name).read_bytes()).hexdigest()
                         for name in sources},
        "method":"Exact U2 initial budget, exact Psi(initial chi), full signed scalar residual profiles, positive residual-weight series.",
        "profiles":[profile_certificate(m) for m in sorted({x[0] for x in CANDIDATES})],
        "finite_candidate_classes":records,
        "assertion_count":len(checks),"checks":checks,
        "limitations":["A failed upper gate is not a proof of divergence or blow-up.",
         "No unbounded search or arbitrary-data regularity claim.",
         "No global Picard convergence follows merely from this endpoint existence gate."]}
Path("enlarged_endpoint_gate_result.json").write_text(json.dumps(result,indent=2)+"\n")
print(result["status"],len(checks),"exact assertions")
