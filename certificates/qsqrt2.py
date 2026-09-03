import os as _os; _ROOT = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))  # release: paths relative to repo root
"""Exact arithmetic in Q(sqrt(2)): elements are (a, b) = a + b*sqrt(2), Fractions."""
from fractions import Fraction as Fr

def zadd(a, b): return (a[0] + b[0], a[1] + b[1])
def zsub(a, b): return (a[0] - b[0], a[1] - b[1])
# (x + y*sqrt2)(u + v*sqrt2) = (xu + 2*yv) + (xv + yu)*sqrt2
def zmul(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def zneg(a): return (-a[0], -a[1])
def zr(a): return (Fr(a), Fr(0))

def zsign(a):
    """Exact sign of a + b*sqrt(2) in {−1,0,+1}."""
    x, y = a
    if x == 0 and y == 0: return 0
    if y == 0: return 1 if x > 0 else -1
    if x == 0: return 1 if y > 0 else -1
    # compare |x| with |y|*sqrt(2) via x^2 vs 2*y^2 (exact Fractions)
    lhs, rhs = x*x, 2*y*y
    if lhs > rhs: return 1 if x > 0 else -1
    if lhs < rhs: return 1 if y > 0 else -1
    raise AssertionError("x^2 == 2 y^2 with nonzero rationals — impossible")

def zgeq(a): return zsign(a) >= 0
def zleq(a): return zsign(a) <= 0

SQRT2 = (Fr(0), Fr(1))
ONE = (Fr(1), Fr(0))
ZERO = (Fr(0), Fr(0))
D4 = (Fr(-1, 4), Fr(1, 4))      # (sqrt2 - 1)/4
D3 = (Fr(-1, 3), Fr(1, 3))      # (sqrt2 - 1)/3

def zrepr(a):
    return f"({a[0]}, {a[1]})*sqrt2-part"


def zsum(iterable, start=None):
    acc = ZERO if start is None else start
    for x in iterable:
        acc = zadd(acc, x)
    return acc


# --- Self-test (runs on import; guards every importer of this module) ---
_S2M1 = (Fr(-1), Fr(1))                    # sqrt2 - 1
assert zmul(SQRT2, SQRT2) == (Fr(2), Fr(0)), "zmul: sqrt2*sqrt2 != 2"
assert zmul(_S2M1, _S2M1) == (Fr(3), Fr(-2)), "zmul: (sqrt2-1)^2 != 3 - 2*sqrt2"
assert zmul(D4, D4) == (Fr(3, 16), Fr(-1, 8)), "zmul: ((sqrt2-1)/4)^2 wrong"
assert zmul(D4, D3) == zmul(D3, D4), "zmul not commutative"
assert zadd(zmul(zr(2), D4), zmul(zr(3), D4)) == zmul(zr(5), D4), "zmul/distributivity"
