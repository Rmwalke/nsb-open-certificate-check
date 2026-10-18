"""One authorized N=3 validation, using new code and exact input coefficients.

All decisions use exact outward integer intervals. Floats occur only in
progress displays. This file and its JSON outputs are intermediate research
records, not a formal or independently certified theorem implementation.
"""
import argparse
import hashlib
import json
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
from n3_interval_core import Context, Solver


class Diagnostics:
    def __init__(self, c, data):
        self.c = c
        self.mass = [F(x) for x in data['norms']]
        self.shells = data['shells']
        self.basis = [{(tuple(v['k']), v['c']): F(v['v'])
                       for v in rows} for rows in data['basis']]
        self.pairs = []
        for row in data['residual_pairs']:
            field = {(tuple(v['k']), v['c']): F(v['v'])
                     for v in row['fourier']}
            norm2 = sum(v*v for v in field.values())
            self.pairs.append((row['i'], row['j'], field,
                               c.sqrt(c.iv(norm2))))
        quartic = defaultdict(F)
        for a, (i,j,f,_) in enumerate(self.pairs):
            for b in range(a, len(self.pairs)):
                k,l,g,_ = self.pairs[b]
                inner = sum(v*g.get(key,0) for key,v in f.items())
                quartic[tuple(sorted((i,j,k,l)))] += inner*(1 if a==b else 2)
        self.quartic = [(tuple(Counter(key).items()), value)
                        for key,value in quartic.items() if value]
        self.origin = []
        self.A1 = []
        for field in self.basis:
            origin = [F(0)]*9
            mode_norm2 = defaultdict(F)
            for (k,comp),v in field.items():
                for deriv in range(3):
                    origin[3*comp+deriv] -= k[deriv]*v
                mode_norm2[k] += v*v
            a1 = c.zero
            for k,norm2 in mode_norm2.items():
                a1 = c.add(a1,c.sqrt(c.iv(sum(t*t for t in k)*norm2)))
            self.A1.append(a1)
            self.origin.append(origin)

    def energy(self, X):
        c=self.c; out=c.zero
        for m,x in zip(self.mass,X):out=c.add(out,c.mulrat(c.sqr(x),m))
        return out

    def residual_square_polynomial(self,X):
        c=self.c
        powers=[]
        for x in X:
            x2=c.sqr(x)
            powers.append([c.one,x,x2,c.mul(x2,x),c.sqr(x2)])
        out=c.zero
        for exponents,coeff in self.quartic:
            value=c.one
            for i,p in exponents:value=c.mul(value,powers[i][p])
            out=c.add(out,c.mulrat(value,coeff))
        if out[1]<0:raise ArithmeticError('negative squared residual upper bound')
        return max(0,out[0]),out[1]

    def residual_square(self,X):
        c=self.c
        direct=self.residual_square_polynomial(X)
        mid=tuple(((x[0]+x[1])//2,)*2 for x in X)
        radius=tuple((max(m[0]-x[0],x[1]-m[0]),)*2 for m,x in zip(mid,X))
        midnorm=c.sqrt(self.residual_square_polynomial(mid))
        error=c.zero
        for i,j,_,norm in self.pairs:
            mi=(abs(mid[i][0]),)*2;mj=(abs(mid[j][0]),)*2
            term=c.add(c.add(c.mul(mi,radius[j]),c.mul(mj,radius[i])),
                       c.mul(radius[i],radius[j]))
            error=c.add(error,c.mul(norm,term))
        lo=max(0,midnorm[0]-error[1]);hi=midnorm[1]+error[1]
        alternative=c.sqr((lo,hi))
        return c.intersect(direct,alternative)

    def gradient(self,X):
        c=self.c
        base=c.sqrt(c.mulrat(c.add(c.sqr(X[0]),c.sqr(X[1])),3))
        upper=base
        for i in range(2,10):
            ax=(max(abs(X[i][0]),abs(X[i][1])),)*2
            upper=c.add(upper,c.mul(ax,self.A1[i]))
        origin2=c.zero
        for q in range(9):
            value=c.zero
            for i in range(10):value=c.add(value,c.mulrat(X[i],self.origin[i][q]))
            origin2=c.add(origin2,c.sqr(value))
        l2=c.zero
        for i in range(10):
            l2=c.add(l2,c.mulrat(c.sqr(X[i]),self.mass[i]*self.shells[i]))
        lower=max(c.sqrt(origin2)[0],c.sqrt(l2)[0])
        if lower>upper[1]:raise ArithmeticError('gradient enclosure conflict')
        return lower,upper[1]


def run(args):
    source=Path('n3_exact_algebra_intermediate.json').read_bytes()
    data=json.loads(source)
    c=Context(args.precision)
    h=F(1,args.denominator)
    solver=Solver(c,data,h=h,order=args.order)
    diag=Diagnostics(c,data)
    X=solver.initial()
    lam=c.zero;eta=c.zero;integ=c.zero;causal=c.zero;Rlower=c.zero
    Cstar=c.iv(27*F(25,3)**6)
    records=[];start=time.monotonic()
    count=int(F(args.horizon)/h)
    result={'kind':'new finite validation; no historical code run',
            'source_sha256':hashlib.sha256(source).hexdigest(),
            'precision_bits':args.precision,'h':str(h),'order':args.order,
            'planned_horizon':args.horizon,'normalization':'(2*pi)^-3 dx',
            'CS':'25/3','nu':1,'cutoff':3,'dimension':10,
            'residual_quartic_terms':len(diag.quartic)}
    def display(x):return [float(c.F(t)) for t in x]
    def record(t,slab,fsq,G,H):
        width=max(x[1]-x[0] for x in X)
        row={'t':str(t),'X':X,'lambda':lam,'eta_prefix':eta,
             'I_prefix':integ,'causal_prefix':causal,'R_lower':Rlower,
             'scalar_obstruction':H,'residual_square':fsq,'gradient':G,
             'max_endpoint_width':str(c.F(width))}
        records.append(row)
        print(json.dumps({'t':str(t),'width':float(c.F(width)),
              'lambda':display(lam),'eta':display(eta),'I':display(integ),
              'H_lower':float(c.F(H[0])),
              'energy':display(diag.energy(X)),
              'elapsed_seconds':round(time.monotonic()-start,2)}),flush=True)
    try:
        for step in range(count):
            t=step*h;T=t+h
            Xnew,slab=solver.step(X)
            fsq=diag.residual_square(slab)
            f=c.mulrat(fsq,2)
            G=diag.gradient(slab)
            growth=c.mulrat(c.add(G,c.sqr(G)),2)
            nextlam=c.add(lam,c.mulrat(growth,h))
            # Λ is nondecreasing. These endpoint extrema enclose the full slab.
            Lex=c.sub((lam[0],nextlam[1]),c.iv(t/2,T/2))
            eta_step=c.mulrat(c.mul(c.exp(c.neg(Lex)),f),h)
            I_step=c.mulrat(c.exp(c.mulrat(Lex,2)),h)
            nexteta=c.add(eta,eta_step)
            # F(s) lies between old and new cumulative forcing, both nonnegative.
            causal=c.add(causal,c.mul(I_step,c.sqr((eta[0],nexteta[1]))))
            eta=nexteta;integ=c.add(integ,I_step)
            # Direct lower recurrence avoids subtracting broad Λ enclosures.
            elo=(growth[0],growth[0])
            first=c.mul(c.exp(c.mulrat(elo,h)),Rlower)
            minrate=min(c.F(growth[0]),F(1,2))
            src=c.mulrat(c.mul(c.exp_rat(t/2+h*minrate),(f[0],f[0])),h)
            Rint=c.add(first,src);Rlower=(Rint[0],Rint[0])
            H=c.mul(c.mul(Cstar,c.exp_rat(-T)),c.sqr(Rlower))
            X=Xnew;lam=nextlam
            if (step+1)%args.progress==0 or step==0 or H[0]>c.Q or step+1==count:
                record(T,slab,fsq,G,H)
                Path('n3_run_progress.json').write_text(json.dumps({'metadata':result,'records':records}))
            if H[0]>c.Q:
                delay=F(1,64)
                while c.mul(H,c.sub(c.one,c.exp_rat(-delay)))[0]<=c.Q:
                    delay*=2
                result.update(status='SCALAR_CONTROL_FAIL',completed_horizon=str(T),
                              scalar_obstruction_lower=str(c.F(H[0])),
                              scalar_blowup_before=str(T+delay),
                              scope='exact H1 scalar majorant with CS=25/3; not PDE failure')
                break
        else:
            T=count*h
            E=diag.energy(X)
            energy_general=c.mulrat(c.exp_rat(-2*T),F(3,2))
            EU=(min(E[1],energy_general[1]),)*2
            # Retained support108: Γ²=660; β*N=18sqrt3.
            Htail=c.sqrt(c.mulrat(EU,660))
            delta=c.add(c.mulrat(Htail,2),c.sqr(Htail))
            Rtail2=c.mulrat(c.sqr(EU),972)
            Lend=c.sub(lam,c.iv(T/2))
            Dtail=c.mulrat(c.mul(Rtail2,c.exp(c.neg(Lend))),F(4,7))
            Itail=c.exp(c.add(c.mulrat(Lend,2),c.mulrat(delta,2)))
            total_eta=c.add(eta,Dtail);total_I=c.add(integ,Itail)
            kappa=c.mul(Cstar,c.mul(c.sqr(total_eta),total_I))
            # β=7/2; exact integral of e^-τ[F+D(1-e^-βτ)]².
            ft=(eta[1],eta[1]);dt=(Dtail[1],Dtail[1])
            bracket=c.add(c.add(c.sqr(ft),c.mulrat(c.mul(ft,dt),F(14,9))),
                          c.mulrat(c.sqr(dt),F(49,72)))
            causal_total=c.mul(Cstar,c.add(causal,c.mul(Itail,bracket)))
            status='PASS_ORIGINAL_GATE' if kappa[1]<c.Q else (
                   'PASS_CAUSAL_GATE' if causal_total[1]<c.Q else 'INCONCLUSIVE')
            result.update(status=status,completed_horizon=str(T),
                 kappa_original_upper=str(c.F(kappa[1])),
                 kappa_causal_upper=str(c.F(causal_total[1])),
                 terminal_energy_upper=str(c.F(EU[1])),
                 eta_tail_upper=str(c.F(Dtail[1])),I_tail_upper=str(c.F(Itail[1])))
    except Exception as error:
        result.update(status='ENCLOSURE_FAILURE',error=repr(error),
                      completed_horizon=str(step*h))
        raise
    finally:
        result['elapsed_seconds']=time.monotonic()-start
        result['records']=records
        Path('n3_validated_result.json').write_text(json.dumps(result,indent=2))
        print('FINAL',json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--precision',type=int,default=224)
    p.add_argument('--denominator',type=int,default=128)
    p.add_argument('--order',type=int,default=16)
    p.add_argument('--horizon',default='8')
    p.add_argument('--progress',type=int,default=32)
    run(p.parse_args())
