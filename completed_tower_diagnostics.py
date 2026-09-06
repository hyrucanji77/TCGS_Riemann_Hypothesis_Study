#!/usr/bin/env python3
"""Complete both xi-moment towers to degree N using Taylor order 4*N+4.

Finite arbitrary-precision diagnostics, NOT interval certificates and NOT a
proof of RH. No zero tables or external datasets enter the computation.
The original diagnostics and shifted_tower_diagnostics scripts are unchanged.
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
    parser.add_argument('--degree', type=int, default=4)
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--scale', type=int, default=100)
    parser.add_argument('--output', type=Path,
                        default=Path('completed_tower_diagnostics.json'))
    args = parser.parse_args()
    if not 0 <= args.degree <= 8 or args.dps < 40 or args.scale <= 0:
        parser.error('Require degree 0..8, dps >= 40, and positive scale.')
    mp.mp.dps = args.dps
    degree = args.degree
    max_moment = 2 * degree + 1
    order = 4 * degree + 4
    a = mp.taylor(xi, mp.mpf('0.5'), order)
    f = [(-1)**j * a[2*j] / a[0] for j in range(max_moment + 2)]
    m = []
    for n in range(max_moment + 1):
        m.append(-(n+1)*f[n+1]
                 - mp.fsum(f[k]*m[n-k] for k in range(1,n+1)))
    scale = mp.mpf(args.scale)
    b = [scale**(n+1)*v for n,v in enumerate(m)]
    def hankel_det(shift: int, N: int):
        return mp.det(mp.matrix([[b[i+j+shift] for j in range(N+1)]
                                 for i in range(N+1)]))
    unshifted = [hankel_det(0,N) for N in range(degree+1)]
    shifted = [hankel_det(1,N) for N in range(degree+1)]
    fmt = lambda x: mp.nstr(x,args.dps-10)
    block_checks = []
    for n in range(1,2*degree+3):
        G = mp.matrix([[2*b[(i+j)//2] if (i+j)%2 == 0 else 0
                        for j in range(n)] for i in range(n)])
        even,odd = (n+1)//2,n//2
        lhs = mp.det(G)
        rhs = 2**n * unshifted[even-1] * (shifted[odd-1] if odd else 1)
        denom = max(abs(lhs),abs(rhs))
        residual = abs(lhs-rhs)/denom if denom else mp.mpf(0)
        if residual > mp.power(10,-(args.dps-15)):
            raise ArithmeticError(f'Block determinant discrepancy at size {n}.')
        block_checks.append({'size':n,'scaled_full_determinant':fmt(lhs),
                             'relative_block_residual':mp.nstr(residual,8)})
    result = {
        'status':'finite floating-point diagnostics; not interval-certified; RH not proved',
        'python':platform.python_version(),'mpmath':mp.__version__,
        'decimal_precision':args.dps,'maximum_tower_degree':degree,
        'maximum_moment_index':max_moment,'taylor_order':order,
        'scale':args.scale,'xi_half':fmt(a[0]),'moments':list(map(fmt,m)),
        'scaled_unshifted_determinants':list(map(fmt,unshifted)),
        'scaled_shifted_determinants':list(map(fmt,shifted)),
        'unshifted_scaling_exponent':'(N+1)^2',
        'shifted_scaling_exponent':'(N+1)(N+2)',
        'block_checks':block_checks,
        'all_computed_prefix_values_positive':all(v>0 for v in unshifted+shifted),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
