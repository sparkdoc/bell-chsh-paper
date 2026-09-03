#!/usr/bin/env python3
"""CGLMP d=3 quantum value: maximally entangled state vs ADGL non-maximally entangled state.
Tries the standard CGLMP measurement phases (Alice alpha=(0,1/2), Bob beta=(1/4,-1/4)) in both
Bob sign conventions and reports the one reproducing CGLMP's 2.87293."""
import numpy as np, math, itertools
T = (12 + 8*math.sqrt(3))/9
def value(gamma, bobsign):
    psi = np.zeros(9, complex)
    for j, c in enumerate([1, gamma, 1]): psi[j*3+j] = c
    psi /= np.linalg.norm(psi)
    def basis(shift, sign):
        return [np.array([np.exp(1j*2*math.pi*j*(sign*k + shift)/3) for j in range(3)])/math.sqrt(3) for k in range(3)]
    A = [basis(0.0, 1), basis(0.5, 1)]; B = [basis(0.25, bobsign), basis(-0.25, bobsign)]
    def P(a, b, ka, kb):
        v = np.kron(A[a][ka], B[b][kb]); return abs(np.vdot(v, psi))**2
    def Peq(a, b, shift): return sum(P(a, b, (kb+shift) % 3, kb) for kb in range(3))  # P(A_a = B_b + shift)
    return (Peq(0, 0, 0) - Peq(0, 0, -1) + Peq(1, 0, -1) - Peq(1, 0, 0) + Peq(1, 1, 0) - Peq(1, 1, -1)
            + Peq(0, 1, 0) - Peq(0, 1, 1))
gam = (math.sqrt(11) - math.sqrt(3))/2
for s in (1, -1):
    v1, v2 = value(1.0, s), value(gam, s)
    print(f"Bob sign {s:+d}: max-entangled I3 = {v1:.6f}   ADGL gamma={gam:.5f}: I3 = {v2:.6f}")
    if abs(v1 - T) < 1e-9:
        print(f"  -> reproduces CGLMP's T = {T:.6f}; ADGL value 1+sqrt(11/3) = {1+math.sqrt(11/3):.6f}; "
              f"F*(N>=4) would be {(v2-2)/8:.6f} instead of {(T-2)/8:.6f}")
        # also scan gamma to confirm the maximum over this family
        gs = np.linspace(0.5, 1.0, 5001); vals = [value(g, s) for g in gs]
        i = int(np.argmax(vals)); print(f"  scan over gamma: max I3 = {vals[i]:.6f} at gamma = {gs[i]:.4f}")
