"""
This file defines the boundaries classes
"""

import scipy as sc
import numpy as np
from math import *

from .utils import coefBinomial



class Boundary:
    """This is a class to represent the boundary condition.
    """
    def __init__(self, B):
        self.B = B

    def __call__(self, r):
        return self.B



class Dirichlet(Boundary):
    """This is a class to represent the homogeneous dirichlet boundary condition.

    :param d: Order of consistency 
    :type d: int
    """
    def __init__(self):
        self.d = 0

    def __call__(self, r, sigma = 0):
        """Builds the left boundary condition

        :param r: Number of left ghost points to define the numerical scheme considered
        :type r: int
        :return: The boundary condition B and a function gn such as (U_{-r},...,U_{-1}) = B (U_0,...,U_{m-1}) + gn(t)
        :rtype: (numpy.ndarray,function)
        """
        return np.zeros((r, 1)), lambda t: 0

    def name(self):
        """Name"""
        return "Dirichlet"


class SILW(Boundary):
    """
    This is a class to represent the Simplified Inverse Lax Wendroff boundary procedure, see [Vilar,Shu, 2015 : Development and stability analysis of the inverse Lax Wendroff boundary treatment for central compact schemes]

    :param kd: Index of the beginning of extrapolation using
    :type kd: int
    :param d: Order of consistency of the boundary condition 
    :type d: int
    :param dx: Space discretization, defaults to 1/10
    :type dx: float, optional
    :param a: Velocity, defaults to 1
    :type a: float, optional
    :param gderivatives: List of the boundary data and its derivatives, defaults to None
    :type gderivatives: list, optional
    """
    
    def __init__(self, kd, d, dx = 0.1, a = 1, gderivatives = None):
        """Constructor method"""
        self.kd = kd
        self.d = d
        self.dx = dx
        self.a = a
        if gderivatives == None:
            gderivatives = [lambda x:0 for i in range(kd)]
        if len(gderivatives) < kd:
            raise IndexError("gderivatives needs at least kd functions")
        self.gderivatives = gderivatives

    def __call__(self, r, sigma = 0):
        """Builds the left boundary condition

        .. note:: m is the number of points needed to define the boundary condition, it is equal to d here.

        :param r: Number of left ghost points to define the numerical scheme considered
        :type r: int
        :param sigma: Gap between the mesh and the boundary condition, defaults to 0
        :type sigma: float, optional
        :return: The boundary condition B and a function gn such as (U_{-r},...,U_{-1}) = B (U_0,...,U_{m-1}) + gn(t)
        :rtype: (numpy.ndarray,function)
        """
        B = np.zeros((r, self.d))
        for j in range(r):
            for m in range(self.kd, self.d):
                for l in range(m + 1):
                    B[r - 1 - j, l] += (-(j + 1) + sigma) ** m / factorial(m) * coefBinomial(m, l) * (-1) ** (m - l)
        def gn_func(t):
            gn = np.zeros(r)
            for j in range(r):
                for k in range(self.kd):
                    gn[r - 1 - j] += (
                        (-(j + 1) + sigma) ** k / factorial(k) * self.dx**k * (-1) ** k * self.gderivatves[k](t) / (self.a**k)
                    )
            return gn
        return B, gn_func

    def name(self):
        """Name
        """
        return f"S{self.kd}ILW{self.d}"


class DDJ(Boundary):
    """
    This is a class to represent the boundary Reconstruction procedure explained in [Boutin, Le Barbenchon, Seguin, 2023 : Stability of finite difference schemes for the hyperbolic initial boundary value problem by winding number computations] and in [Dakin Desprès Jaouen, 2018 : Inverse Lax–Wendroff boundary treatment for compressible Lagrange-remap hydrodynamics on cartesian grids]
    
    .. warning:: Warning d represents the order and kd represente the truncature, the notation is not in the same order than SILW, it is to respect the order of [DakinDespresJaouen18]
    
    :param kd: Index of the beginning of extrapolation using
    :type kd: int
    :param d: Order of consistency of the boundary condition 
    :type d: int
    :param dx: Space discretization, defaults to 1/10
    :type dx: float, optional
    :param a: Velocity, defaults to 1
    :type a: float, optional
    :param gderivatives: List of the boundary data and its derivatives, defaults to None
    :type gderivatives: list, optional
    """
    def __init__(self, d, kd, dx = 0.1, a = 1, gderivatives = None):
        """Constructor method"""
        self.d = d
        self.kd = kd
        self.dx = dx
        self.a = a
        if gderivatives == None:
            gderivatives = [lambda x:0 for i in range(kd)]
        if len(gderivatives) < kd:
            raise IndexError("gderivatives needs at least kd functions") #à vérifier
        self.gderivatives = gderivatives

    def __call__(self, r, sigma = 0):
        """Builds the left boundary condition

        .. note:: m is the number of points needed to define the boundary condition, it is equal to d here.

        :param r: Number of left ghost points to define the numerical scheme considered
        :type r: int
        :param sigma: Gap between the mesh and the boundary condition, defaults to 0
        :type sigma: float, optional
        :return: The boundary condition B and a function gn such as (U_{-r},...,U_{-1}) = B (U_0,...,U_{m-1}) + gn(t)
        :rtype: (numpy.ndarray,function)
        """
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
        #"""write the boundary condition with gn"""
        return ConditionBord[::-1], lambda x: 0

    def name(self):
        """Name"""
        return f"DDJ{self.d},{self.kd}"

