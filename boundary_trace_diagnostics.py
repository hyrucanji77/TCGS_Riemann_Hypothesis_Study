#!/usr/bin/env python3
"""Exact manufactured checks for boundary reduction and resolvent identities.

Run: python boundary_trace_diagnostics.py --output boundary_trace_diagnostics.json
Uses only the Python standard library. It is not a source-geometry solver or
an arithmetic realization of the Riemann xi function.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import platform

Matrix = list[list[F]]


def matrix(rows: list[list[int | F]]) -> Matrix:
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('A nonempty rectangular matrix is required.')
    return [[F(x) for x in row] for row in rows]


def eye(n: int) -> Matrix:
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a: Matrix, b: Matrix, sign: int = 1) -> Matrix:
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError('Matrix shapes do not agree.')
    return [[x + sign*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def mul(a: Matrix, b: Matrix) -> Matrix:
    if len(a[0]) != len(b):
        raise ValueError('Incompatible multiplication shapes.')
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def scale(a: Matrix, c: F) -> Matrix:
    return [[c*x for x in row] for row in a]


def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]


def inverse(a: Matrix) -> Matrix:
    n = len(a)
    if len(a[0]) != n:
        raise ValueError('Only square matrices can be inverted.')
    e = eye(n)
    work = [a[i][:] + e[i] for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if work[i][j]), None)
        if pivot is None:
            raise ValueError('Singular matrix.')
        work[j], work[pivot] = work[pivot], work[j]
        value = work[j][j]
        work[j] = [x/value for x in work[j]]
        for i in range(n):
            if i != j:
                value = work[i][j]
                work[i] = [x-value*y for x, y in zip(work[i], work[j])]
    return [row[n:] for row in work]


def determinant(a: Matrix) -> F:
    n = len(a)
    if len(a[0]) != n:
        raise ValueError('Only square matrices have a determinant.')
    if n == 1:
        return a[0][0]
    return sum(((-1)**j*a[0][j]*determinant(
        [row[:j]+row[j+1:] for row in a[1:]]) for j in range(n)), F(0))


def trace(a: Matrix) -> F:
    return sum((a[i][i] for i in range(len(a))), F(0))


def require(condition: bool, name: str, checks: dict[str, bool]) -> None:
    if not condition:
        raise ArithmeticError('Failed exact check: ' + name)
    checks[name] = True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('boundary_trace_diagnostics.json'))
    args = parser.parse_args()
    checks: dict[str, bool] = {}
    C = matrix([[2, -1, -1], [-1, 2, 0], [-1, 0, 2]])
    B = matrix([[-1, -1]])
    S = add(matrix([[2, 0], [0, 2]]), scale(mul(transpose(B), B), F(1, 2)), -1)
    T = matrix([[F(1, 2), F(1, 2)], [1, 0], [0, 1]])
    K = inverse(S)
    require(S == matrix([[F(3, 2), F(-1, 2)], [F(-1, 2), F(3, 2)]]),
            'schur_response', checks)
    require(mul(mul(transpose(T), C), T) == S, 'minimizing_extension', checks)
    minors = [determinant([row[:j] for row in C[:j]]) for j in (1, 2, 3)]
    require(minors == [F(2), F(3), F(4)], 'positive_leading_minors', checks)
    require(K == matrix([[F(3, 4), F(1, 4)], [F(1, 4), F(3, 4)]]),
            'inverse_normalization', checks)
    power = eye(2)
    moments: list[F] = []
    for n in range(9):
        power = mul(power, K)
        moments.append(trace(power))
    require(all(v == 1 + F(1, 2**(n+1)) for n, v in enumerate(moments)),
            'finite_power_traces', checks)
    hankel = lambda n: [[moments[i+j] for j in range(n+1)] for i in range(n+1)]
    require(determinant(hankel(1)) == F(1, 8) and determinant(hankel(2)) == 0,
            'rank_two_hankel_termination', checks)
    for w in (F(-3), F(-1, 2), F(0), F(1, 3), F(3)):
        P = add(C, scale(eye(3), w), -1)
        interior = P[0][0]
        spectral_S = add([row[1:] for row in P[1:]],
                         scale(mul(transpose([P[0][1:]]), [P[0][1:]]), 1/interior), -1)
        require(determinant(P) == interior*determinant(spectral_S),
                'spectral_factorization_' + str(w), checks)
        require(determinant(P)/4 == (2-w)*((2-w)**2-2)/4,
                'full_characteristic_polynomial_' + str(w), checks)
        require(determinant(add(eye(2), scale(K, w), -1)) == (1-w)*(1-w/2),
                'boundary_characteristic_polynomial_' + str(w), checks)
    Q = matrix([[F(3, 4), 0], [0, F(1, 4)]])
    require(mul(K, Q) != mul(Q, K), 'noncommuting_control', checks)
    diff = add(K, Q, -1)
    entry_trace_bound = sum((abs(x) for row in diff for x in row), F(0))
    for t in (F(1), F(3, 2), F(4)):
        RK = inverse(add(eye(2), scale(K, t*t)))
        RQ = inverse(add(eye(2), scale(Q, t*t)))
        left = add(mul(K, RK), mul(Q, RQ), -1)
        right = mul(mul(RK, diff), RQ)
        require(left == right, 'resolvent_identity_' + str(t), checks)
        require(abs(2*t*trace(left)) <= 2*t*entry_trace_bound,
                'rational_trace_upper_bound_' + str(t), checks)
    result = {
        'status': 'exact rational manufactured algebra checks; not an RH proof',
        'python_version': platform.python_version(),
        'arithmetic': 'fractions.Fraction',
        'source_geometry_realized': False,
        'uses_zeta_zeros': False,
        'checks': checks,
        'number_of_checks': len(checks),
        'schur_response': [[str(x) for x in row] for row in S],
        'inverse': [[str(x) for x in row] for row in K],
        'moments': [str(x) for x in moments],
        'hankel_determinants': [str(determinant(hankel(n))) for n in (0, 1, 2)]
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
