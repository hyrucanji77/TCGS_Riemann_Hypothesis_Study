#!/usr/bin/env python3
"""Morse/Robin asymptotic diagnostics and exact rational coefficient checks.

Floating-point special-function values are not interval certificates.
The exclusion in the manuscript uses the analytic coefficient B_* < 0,
whose coarse rational bound is independently checked below.
Requires mpmath==1.3.0; no zeta-zero routine is used.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import mpmath as mp


def trim(a):
    a = list(map(Q, a))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else Q(0))
                 + (b[i] if i < len(b) else Q(0))
                 for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([v * c for v in a])


def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i+j] += u*v
    return trim(out)


def derivative(a):
    return trim([i*a[i] for i in range(1, len(a))] or [Q(0)])


def compose(a, b):
    out = [Q(0)]
    for c in reversed(a):
        out = add(mul(out, b), [c])
    return out


def atan_bounds(z: Q, terms: int):
    """Alternating-series enclosure, using consecutive rational partial sums."""
    s = sum(((-1)**j * z**(2*j+1) / (2*j+1)
             for j in range(terms)), Q(0))
    following = s + (-1)**terms * z**(2*terms+1) / (2*terms+1)
    return min(s, following), max(s, following)


def exact_checks():
    checks = []
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)
    k = Q(5, 4)
    delta = [k*k/2-Q(3,32), -k/2, Q(1,16)]
    hp = scale(delta, 2)
    check('Robin h polynomial', hp == [Q(11,8), -Q(5,4), Q(1,8)])
    b = add([ -Q(9,128)+k**3/6-k/24, k/4, -Q(1,16)],
            scale(mul(hp,hp),Q(1,8)))
    check('Robin residual quartic', b == [Q(225,512),-Q(60,512),
                                         Q(90,512),-Q(20,512),Q(1,512)])
    c1 = [-k*k/2, k/2, -Q(1,16)]
    c2 = [-k**3/6+k/24, k/4, -Q(1,16)]
    check('First Whittaker x derivative', derivative(c1) == [k/2, -Q(1,8)])
    check('Second Whittaker x derivative', derivative(c2) == [k/4, -Q(1,8)])
    a = derivative(c1)
    check('Riccati first recurrence',
          derivative(c2) == scale(add(mul([Q(0),Q(1)],derivative(a)),a),Q(1,2)))
    b2 = [Q(1,6),-Q(1),Q(1)]
    b3 = [Q(0),Q(1,2),-Q(3,2),Q(1)]
    half_minus_k = [Q(1,2),-Q(1)]
    g1 = scale(add(compose(b2,[Q(1,2)]), scale(compose(b2,half_minus_k),-1)), Q(1,2))
    g2 = scale(add(compose(b3,[Q(1,2)]), scale(compose(b3,half_minus_k),-1)), -Q(1,6))
    check('Gamma ratio first coefficient, polynomial in k', g1 == [Q(0),Q(0),-Q(1,2)])
    check('Gamma ratio second coefficient, polynomial in k', g2 == [Q(0),Q(1,24),Q(0),-Q(1,6)])
    t1 = (compose(b2,[Q(1,4)])[0]-b2[0])/2
    t2 = -(compose(b3,[Q(1,4)])[0]-b3[0])/6
    check('Xi gamma first coefficient', t1 == -Q(3,32))
    check('Xi gamma second coefficient', t2 == -Q(1,128))
    check('Xi completed second coefficient', t2-Q(1,16) == -Q(9,128))
    check('Shift square-root coefficient', 2*Q(1,8) == Q(1,4))
    check('Shift next square-root coefficient', Q(1,8)**2+2*(-Q(1,128)) == 0)
    tan2 = 2*Q(1,5)/(1-Q(1,5)**2)
    tan4 = 2*tan2/(1-tan2*tan2)
    check('Machin angle tangent identity', (tan4-Q(1,239))/(1+tan4*Q(1,239)) == 1)
    a5 = atan_bounds(Q(1,5),8)
    a239 = atan_bounds(Q(1,239),4)
    pi_lo, pi_hi = 16*a5[0]-4*a239[1], 16*a5[1]-4*a239[0]
    lo, hi = Q(157,50), Q(22,7)
    check('Coarse pi enclosure from Machin series', lo < pi_lo < pi_hi < hi)
    h = lambda p: 2*p*p-5*p+Q(11,8)
    check('h monotone on the pi interval',4*lo-5 > 0)
    check('Positive matched h', h(lo)>0)
    check('h upper bound', h(hi) == Q(2123,392))
    upper = Q(13,64)+5*hi/4-lo*lo+h(hi)**2/8
    lower = Q(13,64)+5*lo/4-hi*hi+h(lo)**2/8
    check('Exact negative B upper bound', upper == -Q(1583907247,768320000) and upper<0)
    check('Exact B lower bound', lower == -Q(20700067791,9800000000))
    return {'count':len(checks), 'passed':checks,
            'B_lower_rational':str(lower),'B_upper_rational':str(upper),
            'pi_lower_rational':str(pi_lo),'pi_upper_rational':str(pi_hi)}


def calculate(dps: int):
    if dps < 40:
        raise ValueError('Use at least 40 decimal working digits.')
    mp.mp.dps = dps
    def out(z):
        return mp.nstr(z, min(dps-10,65))
    half = mp.mpf(1)/2
    k = mp.mpf(5)/4
    x = 4*mp.pi
    h = 2*mp.pi**2-5*mp.pi+mp.mpf(11)/8
    B = mp.mpf(13)/64+5*mp.pi/4-mp.pi**2+h*h/8
    def xi(s):
        return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    central = xi(half)
    W0 = mp.whitw(k,0,x)
    Wp0 = mp.diff(lambda xx:mp.whitw(k,0,xx),x)
    a = half-k
    WpU = (half/x-half)*W0-a*mp.exp(-x/2)*mp.sqrt(x)*mp.hyperu(a+1,2,x)
    derivative_error = abs(Wp0-WpU)
    if derivative_error > mp.power(10,-dps+8):
        raise ArithmeticError('Whittaker derivative formulas disagree.')
    Q0 = 2*x*Wp0-(1+h)*W0
    C = -mp.pi**(-mp.mpf(1)/4)*Q0/central
    def Qmu(mu):
        W = mp.whitw(k,mu,x)
        Wp = mp.diff(lambda xx:mp.whitw(k,mu,xx),x)
        return 2*x*Wp-(1+h)*W
    rows=[]
    for mu_int in (30,100,300,1000):
        mu=mp.mpf(mu_int)
        log_ratio = mp.log(xi(2*mu+half)/central)-mp.log(Qmu(mu)/Q0)
        scaled = mu*mu*(log_ratio-mp.log(C))
        c1=-k*k/2+k*x/2-x*x/16
        c2=-k**3/6+k/24+k*x/4-x*x/16
        logL=mp.log(mp.sqrt(x)/(2*mp.sqrt(mp.pi)))+mp.loggamma(mu)-mu*mp.log(x/4)+k*mp.log(mu)
        w_second=mu*mu*(mp.log(mp.whitw(k,mu,x))-logL-c1/mu)
        shift=mp.mpf(4)
        mup=mp.sqrt(mu*mu+shift/4)
        shift_test=8*mu*(mp.log(-Qmu(mup))-mp.log(-Qmu(mu)))/(mp.log(mu)-mp.log(mp.pi))
        rows.append({'mu':mu_int,'mu_squared_log_ratio_over_C':out(scaled),
                     'B_target':out(B),'Whittaker_second_scaled':out(w_second),
                     'Whittaker_c2_target':out(c2),'shift_four_scaled':out(shift_test)})
    return {'working_decimal_digits':dps, 'mpmath_version':mp.__version__,
            'numerical_status':'ordinary arbitrary-precision diagnostics; not interval certificates',
            'exact_algebra':exact_checks(),
            'matched':{'kappa':'5/4','x':out(x),'h':out(h),'B':out(B),
                       'W_at_zero_order':out(W0),'W_x_at_zero_order':out(Wp0),
                       'Q_at_zero_order':out(Q0),'xi_half':out(central),
                       'C_R':out(C), 'C_R_minus_one':out(C-1),
                       'derivative_formula_absolute_difference':out(derivative_error)},
            'asymptotic_checks':rows}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=80)
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    text=json.dumps(calculate(args.dps),indent=2,ensure_ascii=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text,end='')


if __name__=='__main__':
    main()
