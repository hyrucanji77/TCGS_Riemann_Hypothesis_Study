#!/usr/bin/env python3
"""Numerical diagnostics for the calibrated Morse exclusion.

These values check an analytically proved asymptotic expansion. They are not
interval-certified enclosures and are not used as premises in the proof.
No Riemann-zero data are used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import mpmath as mp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps', type=int, default=60)
    parser.add_argument('--output', type=Path, default=Path('morse_diagnostics.json'))
    args = parser.parse_args()
    if args.dps < 25:
        parser.error('--dps must be at least 25')
    mp.mp.dps = args.dps
    kappa = mp.mpf(9) / 4
    x = 4 * mp.pi

    def xi(s):
        return s * (s-1) / 2 * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)

    xi_half = xi(mp.mpf('0.5'))
    W0 = mp.whitw(kappa, 0, x)
    C = 2 * W0 / (xi_half * mp.pi**mp.mpf('0.25'))
    c = -kappa**2/2 + kappa*x/2 - x**2/16
    delta = -mp.mpf(3)/32 - c
    digits = max(20, args.dps - 8)
    fmt = lambda value: mp.nstr(value, digits)
    rows = []
    for n in (30, 100, 300):
        mu = mp.mpf(n)
        s = 2*mu + mp.mpf('0.5')
        W = mp.whitw(kappa, mu, x)
        leading_W = (mp.sqrt(x)/(2*mp.sqrt(mp.pi)) * mp.gamma(mu)
                     * (x/4)**(-mu) * mu**kappa)
        target = xi(s)/xi_half
        determinant = W/W0
        rows.append({
            'mu': n,
            'mu_times_W_relative_remainder': fmt(mu*(W/leading_W-1)),
            'mu_times_normalized_target_ratio_remainder':
                fmt(mu*(target/determinant/C-1)),
        })
    result = {
        'python_version': platform.python_version(),
        'mpmath_version': mp.__version__,
        'working_decimal_precision': args.dps,
        'interval_certified': False,
        'uses_zeta_zeros': False,
        'kappa': fmt(kappa),
        'x': fmt(x),
        'C_kappa': fmt(C),
        'c_kappa': fmt(c),
        'delta_kappa': fmt(delta),
        'rows': rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
