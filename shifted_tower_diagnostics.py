#!/usr/bin/env python3
"""Recompute the shifted Stieltjes tower from central xi derivatives.

These arbitrary-precision results are NOT interval-certified signs or a proof
of RH. No zeta-zero table or routine is used. Existing diagnostics.py is not
modified. Run with --dps 40 and --dps 60 for the manuscript's printed prefix.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import mpmath as mp
from diagnostics import xi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=3,
                        help='Maximum degree of the shifted tower, 0..8.')
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--scale', type=int, default=100)
    parser.add_argument('--output', type=Path, default=Path('shifted_tower_diagnostics.json'))
    args = parser.parse_args()
    if not 0 <= args.degree <= 8 or args.dps < 40 or args.scale <= 0:
        parser.error('Require degree 0..8, dps >= 40, and positive scale.')
    mp.mp.dps = args.dps
    max_moment = 2 * args.degree + 2
    a = mp.taylor(xi, mp.mpf('0.5'), 2 * (max_moment + 1))
    f = [(-1)**j * a[2*j] / a[0] for j in range(max_moment + 2)]
    moments = []
    for n in range(max_moment + 1):
        moments.append(-(n+1)*f[n+1] - mp.fsum(
            f[k]*moments[n-k] for k in range(1, n+1)))
    b = [mp.mpf(args.scale)**(n+1)*v for n, v in enumerate(moments)]
    unshifted = [mp.det(mp.matrix([[b[i+j] for j in range(N+1)]
                                  for i in range(N+1)]))
                 for N in range(args.degree + 2)]
    shifted = [mp.det(mp.matrix([[b[i+j+1] for j in range(N+1)]
                                for i in range(N+1)]))
               for N in range(args.degree + 1)]
    checks = []
    # Scaled lifted moments: 2*b[j] in even degree, 0 in odd degree.
    for n in range(1, 2*args.degree + 3):
        G = mp.matrix([[2*b[(i+j)//2] if (i+j) % 2 == 0 else 0
                        for j in range(n)] for i in range(n)])
        even = (n + 1)//2
        odd = n//2
        rhs = 2**n * unshifted[even-1] * (shifted[odd-1] if odd else 1)
        lhs = mp.det(G)
        residual = abs(lhs-rhs) / max(abs(lhs), abs(rhs))
        if residual > mp.power(10, -(args.dps-15)):
            raise ArithmeticError(f'Block determinant check unstable at size {n}.')
        checks.append({'size': n, 'scaled_full_determinant': mp.nstr(lhs, args.dps-10),
                       'relative_block_residual': mp.nstr(residual, 8)})
    s = lambda value: mp.nstr(value, args.dps-10)
    out = {
        'status': 'finite floating-point diagnostics; not interval-certified; RH not proved',
        'python': platform.python_version(), 'mpmath': mp.__version__,
        'decimal_precision': args.dps, 'maximum_shifted_degree': args.degree,
        'scale': args.scale, 'xi_half': s(a[0]),
        'moments': list(map(s, moments)),
        'scaled_unshifted_determinants': list(map(s, unshifted)),
        'scaled_shifted_determinants': list(map(s, shifted)),
        'unshifted_scaling_exponent': '(N+1)^2',
        'shifted_scaling_exponent': '(N+1)(N+2)',
        'block_checks': checks,
        'all_displayed_prefix_values_positive': all(x > 0 for x in unshifted+shifted)
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
