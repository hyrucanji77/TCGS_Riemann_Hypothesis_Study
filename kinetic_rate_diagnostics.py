#!/usr/bin/env python3
"""Exact rational checks of the kinetic/rate reduction and tower scaling.

The finite examples check algebraic implementation, not an infinite spectral
identity, and are not a proof of RH. The general deduction is in Section 11.7.
This script uses only the Python standard library.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path


def det(a: list[list[Q]]) -> Q:
    """Exact determinant by elimination with row pivoting, including size zero."""
    a = [[Q(x) for x in row] for row in a]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError('Expected a square matrix.')
    ans = Q(1)
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i] != 0), None)
        if pivot is None:
            return Q(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            ans = -ans
        p = a[i][i]
        ans *= p
        for j in range(i+1, n):
            factor = a[j][i]/p
            for k in range(i+1, n):
                a[j][k] -= factor*a[i][k]
            a[j][i] = Q(0)
    return ans


def checks() -> list[str]:
    passed: list[str] = []
    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        passed.append(name)
    for d, beta, c, kappa, h in [
        (Q(1,2), Q(3), Q(5), Q(7,4), Q(-2)),
        (Q(3), Q(2,3), Q(7,3), Q(5,4), Q(11,2)),
        (Q(4), Q(1,2), Q(13), Q(9,4), Q(0))]:
        a = d*d
        b = beta*d
        k_lag = -2*kappa/b
        check(f'kinetic normalization d={d}', a/(d*d) == 1)
        check(f'leading potential coefficient d={d}', b*b/Q(4)*(2*c/b)**2 == c*c)
        check(f'lower potential coefficient d={d}', b*b*k_lag*(2*c/b) == -4*c*kappa)
        check(f'Robin transport d={d}', (d*h)/b == h/beta)
        check(f'counting constant identity d={d}', Q(4)/(b*(2*c/b)) == 2/c)
        if b == 2:
            check(f'calibrated Lagarias sign d={d}', k_lag == -kappa)
    check('first counting coefficient b=2', Q(1, 2*Q(2)) == Q(1,4))
    atoms = [Q(1,2), Q(1,3), Q(1,5), Q(1,7), Q(1,11)]
    m = [sum((z**(n+1) for z in atoms), Q(0)) for n in range(11)]
    scale = Q(100)
    H = lambda n, shift: [[m[i+j+shift] for j in range(n)] for i in range(n)]
    for N in range(4):
        n = N+1
        lhs = det([[scale**(i+j+2)*m[i+j+1] for j in range(n)] for i in range(n)])
        check(f'Stieltjes determinant scaling N={N}', lhs == scale**((N+1)*(N+2))*det(H(n,1)))
    for n in range(1,9):
        G = [[2*m[(i+j)//2] if (i+j)%2 == 0 else Q(0)
              for j in range(n)] for i in range(n)]
        ne, no = (n+1)//2, n//2
        check(f'Grommer block determinant size={n}', det(G) == 2**n*det(H(ne,0))*det(H(no,1)))
    return passed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('kinetic_rate_diagnostics.json'))
    args = parser.parse_args()
    passed = checks()
    out = {'status': 'exact finite algebraic checks, not a proof of RH',
           'exact_checks_passed': len(passed), 'checks': passed}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
