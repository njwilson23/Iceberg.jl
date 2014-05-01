#! /usr/bin/env python
"""
Finds the approximate λ that minimizes the two-phase classical Stefan equation
in 2D. The interface position then is given by X(t) = 2λ sqrt(αl t)
"""

from math import pi, sqrt, exp
from scipy.special import erf, erfc
from scipy.optimize import minimize_scalar

def eval_stefan(λ, **kw):
    """ Objective function to solve transcendental equation for the stefan
    parameter λ """

    Lf = kw.get("Lf", 1.0)
    cpl = kw.get("cpl", 1.0)
    cps = kw.get("cps", 1.0)
    ρl = kw.get("rhol", 1.0)
    ρs = kw.get("rhos", 1.0)
    κl = kw.get("kapl", 2.0)
    κs = kw.get("kaps", 1.0)
    Tl = kw.get("Tl", 1.0)
    Tm = kw.get("Tm", 0.0)
    Ts = kw.get("Ts", -1.0)

    αl = κl / (ρl * cpl)
    αs = κs / (ρs * cps)
    ν = sqrt(αl / αs)
    Stl = cpl * (Tl-Tm) / Lf
    Sts = cps * (Tm-Ts) / Lf

    res = Stl / (exp(λ*λ) * erf(λ)) \
        - Sts / (ν * exp(ν*ν*λ*λ) * erfc(ν*λ)) \
        - sqrt(pi)*λ

    return res**2

# test case: Glauber's salt values
# reference answer: λ = 0.521, with an objective value of 1.07e-9
fsalt = lambda λ: eval_stefan(λ, Tm=32.0, cpl=3.31, cps=1.76, kapl=5.9e-2, 
                                 kaps=2.16e-3, rhos=1460.0, rhol=1460.0,
                                 Lf=251.21, Ts=25.0, Tl=90.0)

fwater = lambda λ: eval_stefan(λ, Tm=0.0, cpl=4.1868, cps=0.5, kapl=5.664e-3,
                                  kaps=2.16e-3, rhos=1.0, rhol=1.0, Lf=333.4)

#print(minimize_scalar(fwater))
print(minimize_scalar(eval_stefan))

