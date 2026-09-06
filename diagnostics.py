#!/usr/bin/env python3
"""Reproduce finite xi-moment diagnostics; NOT an interval-arithmetic proof.

Run: python diagnostics.py --degree 4 --dps 60 --output diagnostics_60.json
Requires mpmath. No zeta-zero ordinates or external datasets are used.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import mpmath as mp


def xi(s):
    """Riemann xi near 1/2 and on the real axis away from 0 and 1."""
    return s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=4)
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--scale', type=int, default=100)
    parser.add_argument('--output', type=Path, default=Path('diagnostics.json'))
    args = parser.parse_args()
    if not 0 <= args.degree <= 12 or args.dps < 40 or args.scale <= 0:
        parser.error('Use degree 0..12, dps >= 40, and a positive scale.')
    mp.mp.dps = args.dps
    max_moment = 2 * args.degree
    a = mp.taylor(xi, mp.mpf('0.5'), 2 * (max_moment + 1))
    f = [(-1)**j * a[2*j] / a[0] for j in range(max_moment + 2)]
    moments = []
    for n in range(max_moment + 1):
        moments.append(-(n+1) * f[n+1] - mp.fsum(
            f[k] * moments[n-k] for k in range(1, n+1)))
    b = [mp.mpf(args.scale)**(n+1) * v for n, v in enumerate(moments)]
    dets = [mp.det(mp.matrix([[b[i+j] for j in range(N+1)]
                             for i in range(N+1)]))
            for N in range(args.degree+1)]
    constant = mp.besselk(0, 2*mp.pi) / (a[0] * (2*mp.pi)**mp.mpf('0.25'))
    ratios = []
    for s in (20, 80, 320):
        s = mp.mpf(s)
        target = xi(s) / a[0]
        candidate = mp.besselk(s/2-mp.mpf('0.25'), 2*mp.pi) / mp.besselk(0, 2*mp.pi)
        normalized = target / candidate / (constant * s**mp.mpf('2.25'))
        ratios.append({'s': int(s), 'normalized_ratio': mp.nstr(normalized, args.dps-10)})
    result = {
        'status': 'finite non-certified numerical diagnostics, not a proof of RH',
        'python': platform.python_version(), 'mpmath': mp.__version__,
        'degree': args.degree, 'decimal_precision': args.dps, 'scale': args.scale,
        'xi_half': mp.nstr(a[0], args.dps-10),
        'moments': [mp.nstr(x, args.dps-10) for x in moments],
        'scaled_hankel_determinants': [mp.nstr(x, args.dps-10) for x in dets],
        'toy_normalized_ratios': ratios,
        'all_computed_leading_determinants_positive': all(x > 0 for x in dets)
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
