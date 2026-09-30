#!/usr/bin/env python3
"""Finite arithmetic-kernel and normalization diagnostics for Appendix E.

Run: python screw_kernel_diagnostics.py --dps 60 --output screw_kernel_60.json
The integer/rational checks are exact. The mpmath evaluations are numerical
regression tests, NOT interval certificates or a proof of global positivity.
No table of zeta zeros, expected output, or source-derived Gram model is read.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import mpmath as mp


def primes_through(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = bytearray(b'\x01') * (limit + 1)
    flags[0:2] = b'\x00\x00'
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b'\x00' * (((limit - p*p) // p) + 1)
    return [p for p in range(2, limit + 1) if flags[p]]


def prime_powers(limit: int) -> list[tuple[int, int]]:
    """Return (integer prime power, underlying prime), each power once."""
    out = []
    for p in primes_through(limit):
        n = p
        while n <= limit:
            out.append((n, p))
            n *= p
    return sorted(out)


def exp_upper(x: F, degree: int = 32) -> F:
    """Exact Taylor/geometric upper bound for exp(x), x >= 0."""
    if x < 0 or x >= degree + 2:
        raise ValueError('Taylor bound requires 0 <= x < degree + 2.')
    term, total = F(1), F(1)
    for k in range(1, degree + 1):
        term *= x / k
        total += term
    first_omitted = term * x / (degree + 1)
    return total + first_omitted / (1 - x / (degree + 2))


def exact_checks() -> dict:
    checks = []
    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append(name)
    for a, cap in [(F(1, 2), 3), (F(1), 8), (F(2), 55)]:
        check(f'complete window cutoff a={a}, P={cap}', exp_upper(2*a) < cap)
    pairs = prime_powers(55)
    check('unique prime powers', len(set(n for n, _ in pairs)) == len(pairs))
    check('composite non-power absent', all(n != 6 for n, _ in pairs))
    check('power 49 retained', (49, 7) in pairs)
    check('power 32 retained', (32, 2) in pairs)
    for n, p in pairs:
        q = n
        while q % p == 0:
            q //= p
        check(f'integer prime-power factorization {n}', q == 1)
    lam = [F(1, 2), F(1, 3)]
    for n in range(7):
        m = sum(x**(n+1) for x in lam)
        doubled = sum(x**(n+1) for x in lam + lam)
        scaled = sum((x/2)**(n+1) for x in lam + lam)
        check(f'signed-spectrum doubled trace {n}', doubled == 2*m)
        check(f'scalar division has wrong higher trace {n}', scaled == F(1, 2)**n*m)
    # The exact x=36 aggregation regression concerns the 2020 printed ordering;
    # the author's later formula already uses the candidate-first ordering.
    array = [[1, 0, 0, 1], [0, 0, 0, 0]]
    old = sum(int(sum(row) > 0) for row in array)
    later = sum(int(sum(array[i][j] for i in range(2)) > 0) for j in range(4))
    check('2020 row aggregation gives one', old == 1)
    check('later candidate aggregation gives two', later == 2)
    check('exact prime count below 36 is eleven', len(primes_through(35)) == 11)
    check('candidate-first x=36 reconstruction', 11 + 2 - later == 11)
    # Perturb g by a bounded continuous even ramp; g(0) remains fixed exactly.
    points = [F(-1), F(-1, 2), F(1, 2), F(1)]
    delta = F(1, 7)
    def perturb(t: F) -> F:
        return delta * min(abs(t), F(2)) / 2
    error = [[perturb(t-u)-perturb(t)-perturb(-u) for u in points] for t in points]
    check('anchored three-term maximum entry error', max(abs(x) for row in error for x in row) <= 3*delta)
    for k, c in enumerate([[F(1)]*4, [F(1),F(-2),F(3),F(-4)],
                           [F(1,3),F(2,5),F(-1,2),F(1,7)]]):
        q = sum(c[i]*error[i][j]*c[j] for i in range(4) for j in range(4))
        check(f'entrywise-to-quadratic bound {k}', abs(q) <= 3*delta*sum(abs(x) for x in c)**2)
        check(f'Euclidean coefficient bound {k}', abs(q) <= 12*delta*sum(x*x for x in c))
    return {'count': len(checks), 'all_pass': True, 'checks': checks}


def xi(s):
    if s == 0 or s == 1:
        return mp.mpf('0.5')
    return s*(s-1)*mp.power(mp.pi, -s/2)*mp.gamma(s/2)*mp.zeta(s)/2


def run(dps: int) -> dict:
    mp.mp.dps = dps
    text = lambda x: mp.nstr(x, dps - 5)
    pairs = prime_powers(55)
    def psi_ar(u):
        u = abs(mp.mpf(u))
        if not u:
            return mp.mpf(0)
        terms = [mp.log(p)/mp.sqrt(n)*(u-mp.log(n))
                 for n, p in pairs if mp.log(n) <= u]
        return (4*(mp.exp(u/2)+mp.exp(-u/2)-2) - mp.fsum(terms)
                + u*(mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi))/2
                + (mp.zeta(2, mp.mpf(1)/4)
                   - mp.exp(-u/2)*mp.lerchphi(mp.exp(-2*u),2,mp.mpf(1)/4))/4)
    def kernel(t, u):
        return -psi_ar(t-u)+psi_ar(t)+psi_ar(u)
    windows = []
    for a, cap in [(mp.mpf('0.5'),3),(mp.mpf(1),8),(mp.mpf(2),55)]:
        pts = [-a,-a/2,a/2,a]
        # Cache repeated differences: special-function calls dominate runtime.
        values = {text(abs(x)): psi_ar(abs(x)) for x in set(pts+[t-u for t in pts for u in pts])}
        val = lambda x: values[text(abs(x))]
        matrix = mp.matrix([[-val(t-u)+val(t)+val(u) for u in pts] for t in pts])
        eigenvalues = mp.eigsy(matrix,eigvals_only=True)
        windows.append({'a':text(a),'complete_integer_cutoff':cap,
                        'points':[text(x) for x in pts],
                        'matrix':[[text(matrix[i,j]) for j in range(4)] for i in range(4)],
                        'eigenvalues':[text(x) for x in eigenvalues],
                        'finite_matrix_positive_at_working_precision':bool(min(eigenvalues)>0)})
    xi0=xi(mp.mpf(1)/2)
    def entire_F(w):
        return xi(mp.mpf(1)/2+mp.j*mp.sqrt(w))/xi0
    def moment(w):
        return -mp.diff(entire_F,w)/entire_F(w)
    conversions=[]
    for z in [mp.mpc(mp.mpf(1)/3,mp.mpf(1)/4),mp.mpc(2,1),mp.mpc(mp.mpf(3)/2,0)]:
        s=mp.mpf(1)/2-mp.j*z
        ratio=mp.diff(xi,s)/xi(s)
        M=moment(z*z)
        q_error=abs(mp.j*ratio-2*z*M)
        R_error=abs(1/(1+ratio)-1/(1-2*mp.j*z*M))
        if max(q_error,R_error)>mp.power(10,-dps+8):
            raise AssertionError('Meromorphic normalization regression failed.')
        conversions.append({'z':text(z),'Q_minus_2zM_abs':text(q_error),
                            'ratio_identity_abs_error':text(R_error)})
    m0=mp.diff(xi,mp.mpf(1)/2,2)/(2*xi0)
    mu0=1+mp.euler/2-mp.log(4*mp.pi)/2
    shifted=moment(-mp.mpf(1)/4)
    if abs(mu0-shifted)>mp.power(10,-dps+8):
        raise AssertionError('Quarter-shift regression failed.')
    return {'status':'PASS','dps':dps,'python':platform.python_version(),
            'mpmath':mp.__version__, 'exact_checks':exact_checks(),
            'arithmetic_windows':windows,'meromorphic_checks':conversions,
            'central_values':{'m0':text(m0),'mu0':text(mu0),
                              'M_minus_quarter':text(shifted),'m0_minus_mu0':text(m0-mu0)},
            'certificate_scope':'Exact rational checks plus arbitrary-precision regression; not interval certification, global kernel positivity, or source identification.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=60)
    parser.add_argument('--output',type=Path,default=Path('screw_kernel_60.json'))
    args=parser.parse_args()
    if not 30<=args.dps<=200:
        parser.error('dps must be between 30 and 200.')
    result=run(args.dps)
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'dps':args.dps,
                      'exact_checks':result['exact_checks']['count'],
                      'windows':len(result['arithmetic_windows'])}))

if __name__=='__main__':
    main()
