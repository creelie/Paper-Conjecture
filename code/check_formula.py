"""Compare the degrees in data/ with the closed formula, and check the gcd theorem.

Usage: python3 code/check_formula.py   (needs sympy)
Each file data/fPQ.txt holds lines "m deg X_m" computed by run_degrees.sh.
"""
import glob
import os
import sys
from fractions import Fraction
from math import gcd

import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def degree(p, q, m):
    """deg X_m(F_{p,q}) by the main theorem."""
    c = p * q - 1
    r, s = m % (p + 1), m % (q + 1)
    main = Fraction((p * (q + 1) + m * c) * (q * (p + 1) + m * c), (p + 1) * (q + 1))
    return main - c * (Fraction(r * (p - r), p + 1) + Fraction(s * (q - s), q + 1))


def check_data():
    ok = True
    for fn in sorted(glob.glob(os.path.join(ROOT, "data", "f??.txt"))):
        p, q = int(fn[-6]), int(fn[-5])
        rows = [ln.split() for ln in open(fn) if ln.strip()]
        bad = [(int(m), int(v)) for m, v in rows if degree(p, q, int(m)) != int(v)]
        note = "" if p >= 2 and q >= 2 else " (p or q = 1: outside the theorem)"
        print(f"(p,q)=({p},{q}): {len(rows)} values, {'all agree' if not bad else 'MISMATCH ' + str(bad[:3])}{note}")
        ok &= not bad
    return ok


def check_gcd(pmax=4, qmax=4, lmax=10):
    a, b = sp.symbols("a b")
    ok = True
    for p in range(1, pmax + 1):
        for q in range(1, qmax + 1):
            A, B = a, b
            for l in range(1, lmax + 1):
                A = sp.expand(b**p * sp.diff(A, a) + a**q * sp.diff(A, b))
                B = sp.expand(b**p * sp.diff(B, a) + a**q * sp.diff(B, b))
                r, s = (l - 1) % (p + 1), (l - 1) % (q + 1)
                want = a**((q - s) * (s >= 1)) * b**((p - r) * (r >= 1))
                ratio = sp.cancel(sp.gcd(A, B) / want)  # gcd over Z carries an integer content
                if not ratio.is_Number:
                    print("gcd mismatch", p, q, l)
                    ok = False
    print(f"gcd theorem checked for 1<=p<={pmax}, 1<=q<={qmax}, 1<=l<={lmax}: {'ok' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    print("F_3 values:", [int(degree(2, 3, m)) for m in range(11)])
    good = check_data() & check_gcd()
    sys.exit(0 if good else 1)
