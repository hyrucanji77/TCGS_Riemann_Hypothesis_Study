#!/usr/bin/env python3
"""Check the free-endpoint comparison and quarter-shift normalization.

Run: python scale_endpoint_diagnostics.py --dps 60 --output scale_endpoint_diagnostics.json
No Riemann-zero ordinates are used. Special-function values are floating-point
checks, not interval certificates and not premises in the analytical proofs.
The polynomial coefficients and finite Grommer block checks use exact fractions.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import platform
import mpmath as mp


def determinant(a: list[list[Q]]) -> Q:
    """Exact determinant by fraction-preserving Gaussian elimination."""
    if not a:
        return Q(1)
    b = [row[:] for row in a]
    n = len(b)
    if any(len(row) != n for row in b):
        raise ValueError('A square matrix is required.')
    result = Q(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if b[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            b[j], b[pivot] = b[pivot], b[j]
            result = -result
        p = b[j][j]
        result *= p
        for i in range(j + 1, n):
            ratio = b[i][j] / p
            for k in range(j + 1, n):
                b[i][k] -= ratio * b[j][k]
            b[i][j] = Q(0)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--output', type=Path, default=Path('scale_endpoint_diagnostics.json'))
    args = parser.parse_args()
    if args.dps < 40:
        parser.error('--dps must be at least 40')
    mp.mp.dps = args.dps
    fmt = lambda a: mp.nstr(a, args.dps - 10)
    half = mp.mpf('0.5')
    kappa = mp.mpf(9) / 4

    def xi(s):
        return s*(s-1)/2 * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)

    # No evaluation of the uncancelled formula exactly at its removable s=1 point.
    m0 = mp.diff(xi, half, 2)/(2*xi(half))
    shift = 1 + mp.euler/2 - mp.log(4*mp.pi)/2
    rows = []
    for name, x in (('lower_cancellation_root', 9-mp.sqrt(42)),
                    ('counting_calibration', 4*mp.pi),
                    ('upper_cancellation_root', 9+mp.sqrt(42))):
        delta = kappa**2/2-kappa*x/2+x*x/16-mp.mpf(3)/32
        tests = []
        for mu_int in (30, 100, 300):
            mu = mp.mpf(mu_int)
            s = 2*mu + half
            # W(kappa,0,x) and xi(1/2) cancel from this relative expression.
            normalized = (xi(s)*mp.sqrt(x) /
                (4*mp.pi**mp.mpf('.25')*mp.whitw(kappa,mu,x)
                 * (x/(4*mp.pi))**mu * mu**(mp.mpf(9)/4-kappa)))
            tests.append({'mu':mu_int,
                          'mu_times_normalized_ratio_minus_one':fmt(mu*(normalized-1))})
        rows.append({'endpoint':name, 'x':fmt(x), 'delta':fmt(delta),
                     'exponential_mismatch_rate':fmt(mp.log(x/(4*mp.pi))),
                     'counting_second_coefficient_difference':fmt(-mp.log(x/(4*mp.pi))/(2*mp.pi)),
                     'large_order_checks':tests})

    # Exact symbolic coefficients of delta_9/4(x), in ascending powers of x.
    coefficients = [Q(9,4)**2/2-Q(3,32), -Q(9,4)/2, Q(1,16)]
    assert coefficients == [Q(39,16), Q(-18,16), Q(1,16)]
    assert Q(81)-Q(39) == Q(42)

    # Manufactured positive reciprocal nodes, not Riemann zeros.
    nodes = (Q(1,4), Q(1,9), Q(1,25), Q(1,49))
    m = lambda n: sum((x**(n+1) for x in nodes), Q(0))
    g = lambda n: 2*m(n//2) if n%2 == 0 else Q(0)
    block_checks = []
    for n in range(1,9):
        even = (n+1)//2
        odd = n//2
        G = [[g(i+j) for j in range(n)] for i in range(n)]
        H = [[m(i+j) for j in range(even)] for i in range(even)]
        Hp = [[m(i+j+1) for j in range(odd)] for i in range(odd)]
        left = determinant(G)
        right = 2**n * determinant(H) * determinant(Hp)
        assert left == right
        block_checks.append({'size':n, 'exact_determinant_identity':True})
    result = {'python':platform.python_version(), 'mpmath':mp.__version__,
              'working_decimal_precision':args.dps, 'interval_certified':False,
              'uses_zeta_zeros':False,
              'm0':fmt(m0), 'M_minus_quarter':fmt(shift),
              'quarter_shift_difference':fmt(m0-shift),
              'morse_comparisons':rows,
              'exact_delta_polynomial_coefficients':[str(x) for x in coefficients],
              'exact_delta_polynomial_check':True,
              'exact_root_polynomial_check':True,
              'grommer_manufactured_checks':block_checks,
              'total_additional_exact_checks':10}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
