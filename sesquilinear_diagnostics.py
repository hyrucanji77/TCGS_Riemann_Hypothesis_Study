#!/usr/bin/env python3
"""Exact complex-coefficient checks of the manuscript's moment-form convention.

All scalar arithmetic uses pairs of fractions (Gaussian rationals). These
manufactured checks detect a missing conjugation; they do not construct the
source-to-arithmetic map or prove the Riemann hypothesis.
"""
from __future__ import annotations
import argparse
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

@dataclass(frozen=True)
class G:
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other: G) -> G:
        return G(self.re+other.re,self.im+other.im)

    def __neg__(self) -> G:
        return G(-self.re,-self.im)

    def __sub__(self, other: G) -> G:
        return self+(-other)

    def __mul__(self, other: G) -> G:
        return G(self.re*other.re-self.im*other.im,
                 self.re*other.im+self.im*other.re)

    def conj(self) -> G:
        return G(self.re,-self.im)

    def json(self) -> dict:
        return {'real':str(self.re),'imaginary':str(self.im)}

ZERO=G()
ONE=G(F(1))
I=G(F(0),F(1))
LAMBDAS=(F(1),F(1,2))

def total(xs) -> G:
    ans=ZERO
    for x in xs:
        ans=ans+x
    return ans

def moment(n: int) -> F:
    return sum((v**(n+1) for v in LAMBDAS),F(0))

def form(p: tuple[G,...],q: tuple[G,...],bar: bool=True) -> G:
    return total((pi.conj() if bar else pi)*qj*G(moment(i+j))
                 for i,pi in enumerate(p) for j,qj in enumerate(q))

def poly(p: tuple[G,...],t: F) -> G:
    return total(c*G(t**i) for i,c in enumerate(p))

def trace_form(p: tuple[G,...],q: tuple[G,...]) -> G:
    return total(G(t)*poly(p,t).conj()*poly(q,t) for t in LAMBDAS)

def multiply(z: G,p: tuple[G,...]) -> tuple[G,...]:
    return tuple(z*c for c in p)

def add(p: tuple[G,...],q: tuple[G,...]) -> tuple[G,...]:
    return tuple((p[j] if j<len(p) else ZERO)+(q[j] if j<len(q) else ZERO)
                 for j in range(max(len(p),len(q))))

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('sesquilinear_diagnostics.json'))
    args=parser.parse_args()
    p=(G(F(1),F(1)),G(F(2),F(-1)),G(F(-1),F(2)))
    q=(G(F(2),F(-1)),G(F(3),F(2)))
    r=(G(F(-1),F(1)),G(F(1,2),F(-3,2)))
    alpha=G(F(2,3),F(-5,4))
    checks=[]
    def check(name: str,truth: bool) -> None:
        if not truth:
            raise ArithmeticError(name)
        checks.append({'name':name,'passed':True})
    check('nontrivial complex cross form',form(p,q).im!=0)
    check('conjugate symmetry',form(q,p)==form(p,q).conj())
    check('first slot conjugates i',form(multiply(I,p),q)==(-I)*form(p,q))
    check('second slot preserves i',form(p,multiply(I,q))==I*form(p,q))
    check('first slot conjugates general scalar',form(multiply(alpha,p),q)==alpha.conj()*form(p,q))
    check('second slot preserves general scalar',form(p,multiply(alpha,q))==alpha*form(p,q))
    check('additivity first slot',form(add(p,r),q)==form(p,q)+form(r,q))
    check('additivity second slot',form(p,add(q,r))==form(p,q)+form(p,r))
    check('cross trace equals Hermitian moment form',trace_form(p,q)==form(p,q))
    check('trace phase test',trace_form(multiply(I,p),q)==(-I)*trace_form(p,q))
    for j,c in enumerate((p,q,r),1):
        v=form(c,c)
        real=tuple(G(x.re) for x in c)
        imag=tuple(G(x.im) for x in c)
        check(f'{j}: real nonnegative quadratic value',v.im==0 and v.re>=0)
        check(f'{j}: phase invariance',form(multiply(I,c),multiply(I,c))==v)
        check(f'{j}: norm-square scalar law',form(multiply(alpha,c),multiply(alpha,c))==alpha.conj()*alpha*v)
        check(f'{j}: real-imaginary decomposition',v==form(real,real)+form(imag,imag))
        check(f'{j}: trace diagonal equals quadratic form',trace_form(c,c)==v)
        for t in LAMBDAS:
            # Coefficient conjugation in p(K)^* for the positive diagonal K.
            check(f'{j}: adjoint at eigenvalue {t}',poly(tuple(x.conj() for x in c),t)==poly(c,t).conj())
    check('negative control: unbarred first slot is linear',
          form(multiply(I,p),q,False)==I*form(p,q,False))
    check('negative control: unbarred and Hermitian cross forms differ',
          form(p,q,False)!=form(p,q))
    check('negative control: p=i gives negative bilinear square',
          form((I,),(I,),False)==G(-moment(0)))
    check('p=i gives positive Hermitian square',form((I,),(I,))==G(moment(0)))
    result={
        'status':'exact manufactured Gaussian-rational checks; not a proof of RH',
        'inner_product_convention':'conjugate-linear in the first slot',
        'eigenvalues':[str(t) for t in LAMBDAS],
        'cross_form':form(p,q).json(),'quadratic_form':form(p,p).json(),
        'bilinear_i_negative_control':form((I,),(I,),False).json(),
        'hermitian_i_control':form((I,),(I,)).json(),
        'check_count':len(checks),'all_passed':True,'checks':checks}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
