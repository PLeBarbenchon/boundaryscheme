"""
This file defines the boundaries classes
"""

import scipy as sc
import numpy as np
from math import *


class Bord:
    def __init__(self, B):
        self.boundary = B


class Dirichlet(Bord):
    def __init__(self):
        self.d = 0

    def __call__(self, r, **kwargs):
        return np.zeros((r, 1)), lambda t: 0

    def name(self):
        return "Dirichlet"


class SILW(Bord):
    def __init__(self, kd, d):
        self.kd = kd
        self.d = d

    def __call__(self, r, **kwargs):
        """
        d = m
        return a function that builds the matrix
        b_{-r,0} ... b_{-r,m-1}
           |            |
        b_{-1,0} ... b_{-1,m-1}
        and a function gn such as (U_{-r},...,U_{-1}) = B (U_0,...,U_{m-1}) + gn(t)
        """
        sigma = kwargs.get("sigma", 0)
        dx = kwargs.get("dx", 0.1)
        a = kwargs.get("a", 1)
        g = kwargs.get("g", lambda x: 0)
        dg = kwargs.get("dg", lambda x: 0)
        d2g = kwargs.get("d2g", lambda x: 0)
        d3g = kwargs.get("d3g", lambda x: 0)
        B = np.zeros((r, self.d))
        for j in range(r):
            for m in range(self.kd, self.d):
                for l in range(m + 1):
                    B[r - 1 - j, l] += (-(j + 1) + sigma) ** m / factorial(m) * comb(m, l) * (-1) ** (m - l)
        functions = [g, dg, d2g, d3g]

        def gn_func(t):
            gn = np.zeros(r)
            for j in range(r):
                for k in range(self.kd):
                    gn[r - 1 - j] += (
                        (-(j + 1) + sigma) ** k / factorial(k) * dx**k * (-1) ** k * functions[k](t) / (a**k)
                    )
            return gn

        return B, gn_func

    def name(self):
        return f"S{self.kd}ILW{self.d}"


class DDJ(Bord):
    """
    Warning d represents the order  and kd represente the truncature, the notation is not in the same order than SILW, it is to respect the order of [DakinDespresJaouen18]
    """

    def __init__(self, d, kd):
        self.d = d
        self.kd = kd

    def __call__(self, r, **kwargs):
        """
        renvoie une fonction qui construit la matrice
        return a function that build the matrix
        b_{-r,0} ... b_{-r,m-1}
           |            |
        b_{-1,0} ... b_{-1,m-1}
        and the function gn
        """
        sigma = kwargs.get("sigma", 0)
        dx = kwargs.get("dx", 0.1)
        a = kwargs.get("a", 1)
        g = kwargs.get("g", lambda x: 0)
        dg = kwargs.get("dg", lambda x: 0)
        d2g = kwargs.get("d2g", lambda x: 0)
        d3g = kwargs.get("d3g", lambda x: 0)
        ymoins = np.zeros((r, self.d - self.kd - 1))
        yplus = np.zeros((self.d - self.kd - 1, self.d - self.kd - 1))
        for i in range(1, r + 1):
            for j in range(1, self.d - self.kd):
                ymoins[i - 1, j - 1] = (
                    (1 - i + 1 / 2 - sigma) ** (self.kd + j + 1) - (1 - i - 1 / 2 - sigma) ** (self.kd + j + 1)
                ) / (factorial(self.kd + j + 1))
        for i in range(1, self.d - self.kd):
            for j in range(1, self.d - self.kd):
                yplus[i - 1, j - 1] = (
                    (i + 1 / 2 - sigma) ** (self.kd + j + 1) - (i - 1 / 2 - sigma) ** (self.kd + j + 1)
                ) / (factorial(self.kd + j + 1))
        if np.linalg.det(yplus) == 0:
            print(yplus, sigma)
        ConditionBord = ymoins.dot(np.linalg.inv(yplus))
        """Ecrire la condition de bord avec le gn"""
        return ConditionBord[::-1], lambda x: 0

    def name(self):
        return f"DDJ{self.d},{self.kd}"

