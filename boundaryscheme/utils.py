"""
This file defines utils function specific to the library
"""

import numpy as np
from math import *
from cmath import *


def coefBinomial(n, k):
    """Computes the binomial coefficient

    :param n: Number of elements
    :type n: int
    :param k: Number of choosen elements
    :type k: int
    :return: Binomial coefficient "n choose k"
    :rtype: int
    """
    if k > (n - k):
        k = n - k
    coef = 1
    for i in range(k):
        coef = coef * (n - i)
        coef = coef // (i + 1)
    return coef


def sort(L):
    """Sort (bubble sort)

    :param L: List of complex numbers
    :type L: list
    :return: The same list but in the ascending order of modulus and for values with same modulus, sort in the ascending order of argument (between -pi and pi)
    :rtype: list
    """
    nbr = len(L)
    for i in range(nbr):
        for j in range(nbr - 1):
            if abs(L[j + 1]) < abs(L[j]):
                x = L[j + 1]
                L[j + 1] = L[j]
                L[j] = x
    for i in range(nbr):
        for j in range(nbr - 1):
            if abs(L[j + 1]) == abs(L[j]):
                if phase(L[j + 1]) < phase(L[j]):
                    x = L[j + 1]
                    L[j + 1] = L[j]
                    L[j] = x
    return L


def sector(z):
    """Finds the sector of z

    .. warning:: the argument is between -pi and pi and here "a" is between 0 and 7

    :param z: Complex number
    :type z: complex
    :return: The sector of z, it is an integer "a" between 0 and 7 which correspond to the angular sector [a*pi/4, (a+1)*pi/4[
    :rtype: int
    """
    ph = phase(z) * 4 / pi
    if ph < 0:
        return int(ph) + 7
    else:
        return int(ph)


def neighboor(sector1, z2):
    """Checks if z2 is in a neighboring sector of sector1

    :param sector1: Integer between 0 and 7 and represent the angular sector [sector1*pi/4, (sector1+1)*pi/4[
    :type sector1: int
    :param z2: Complex number
    :type z2: complex
    :return: True if and only if the complex z2 is in a neighboring sector of [sector1*pi/4, (sector1+1)*pi/4[
    :rtype: bool
    """
    return abs((sector(z2) - sector1) % 8) <= 1


def parametrization(n_param, curve_formula):
    """Finds the discretization of the curve refine if it is needed (as [ZapataMartin2014] procedure).

    :param n_param: Integer for the default discretization of [0,2pi]
    :type n_param: int
    :param curve_formula: A fonction : z in the unit circle mapsto a complex and represent a curve
    :type curve_fomula: function
    :return: The discretization of the unit circle and the values of the curve for this discretization
    :rtype: (list,list)
    """
    dx = 2 * pi / n_param
    current_dx = dx
    current_param = 0
    Param = [current_param]
    curve = [curve_formula(1)]
    current_sector = sector(curve_formula(np.exp(1j * Param[0])))  # calcul le secteur du premier point
    c = 0
    s = 0
    while Param[-1] < 2 * pi:
        current_param = Param[-1] + current_dx
        current_curve_point = curve_formula(np.exp(1j * current_param))
        if neighboor(current_sector, current_curve_point):
            Param.append(current_param)
            curve.append(current_curve_point)
            current_sector = sector(current_curve_point)
            current_dx = dx
            c = 0
        elif c < 40:
            current_dx = current_dx / 2
            c += 1
        else:
            current_param = Param[-1] + dx
            Param.append(current_param)
            curve.append(curve_formula(np.exp(1j * current_param)))
            current_sector = sector(curve[-1])
            current_dx = dx
            c = 0
        s += 1
    if Param[-1] > 2 * pi:
        Param[-1] = 0
        curve[-1] = curve_formula(1)
    return np.array(Param), np.array(curve)


def epsilon(L):
    """Computes the half of the smallest distance between two distincts elements of the list L

    :param L: List of complex numbers (z_i)
    :type L: list
    :return: min \|z_i-z_j| / 2 (when z_i != z_j)
    :rtype: float
    """
    assert len(L) > 1
    diff = []
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if L[i] != L[j]:
                diff.append(abs(L[i] - L[j]))
    mini = min(diff)
    assert mini > 0
    return mini / 2


def eta_func(eps, kappa0, N, polynom, r):
    """Computes the distance eta
    .. note:: for more details, see [Boutin, Le Barbenchon, Seguin, 2023 : Stability of finite difference schemes for the hyperbolic initial boundary value problem by winding number computations]

    :param eps: Radius
    :type eps: float
    :param kappa0: Center
    :type kappa0: complex
    :param N: Number of discretization of the circle centered in kappa0 of radius eps
    :type N: int
    :param polynom: Polynomial
    :type polynom: Polynomial
    :param r: Integer
    :type r: int
    :return: min \|polynom(kappa)|/(1+eps)^r for kappa on the circle centered in kappa0 of radius eps
    :rtype: float
    """
    theta = np.linspace(0, 2 * pi, N)
    circle = np.cos(theta) + 1j * np.sin(theta)
    kappas = kappa0 + eps * circle
    val_pol = np.abs(polynom(kappas))
    mini = min(val_pol)
    assert mini > 0
    return mini / ((1 + eps) ** r)


def schemsum(schem1, schem2):
    """Addition of two numerical schemes

    :param schem1: Representation of a numerical scheme with the interior coefficients and the center
    :type schem1: (list,int)
    :param schem2: Representation of a numerical scheme with the interior coefficients and the center
    :type schem2: (list,int)
    :return: Sum of schem1 and schem2
    :rtype: (list,int)
    """
    inter1, center1 = schem1
    inter2, center2 = schem2
    center = max(center1, center2)
    p = max(len(inter1) - center1, len(inter2) - center2) - 1
    newinter = np.zeros(p + center + 1)
    for i in range(len(inter1)):
        newinter[center - center1 + i] += inter1[i]
    for i in range(len(inter2)):
        newinter[center - center2 + i] += inter2[i]
    return newinter, center


def schemscal(schem, scal):
    """Multiplication of a scheme by a scalar

    :param schem: Representation of a numerical scheme with the interior coefficients and the center
    :type schem: (list,int)
    :param scal: Scalar element
    :type scal: float
    :return: Multiplication of schem by scal
    :rtype: (list,int)
    """
    inter, center = schem
    return scal * inter, center


def boundary_to_boundary_quasi_toep(boundary_condition, gn, inter, center):
    """Transformation from a boundary condition written with ghost points to a boundary condition written as the extraction of the r first rows of the Quasi-Toeplitz matrix

    :param boundary_condition: Boundary condition written as U_{-r} = ..., ...,U_{-1} =..., etc
    :type boundary_condition: numpy.ndarray
    :param gn: Boundary data
    :type gn: function
    :param inter: List of float coefficient that represent the interior scheme
    :type inter: list
    :param center: Index of the central coefficient of the scheme
    :type center: int
    :return: Extraction of the r first rows of the Quasi-Toeplitz matrix and the function b_n(t) such that (U_0^{n+1},...,U_{r-1}^{n+1}) = T_J (U_0^n,...,U_{m-1}^n) + b_n(t^n)
    :rtype: (numpy.ndarray,function)
    """
    r = center
    p = len(inter) - r - 1
    m = max(len(boundary_condition[0]), r + p)
    A = np.zeros((r, r))
    for j in range(r):
        A += inter[j] * np.diag(np.ones(r - j), j)
    Ap = np.zeros((r, m))
    for j in range(1, r + p + 1):
        if j < r:
            Ap[:r, :r] += inter[j] * np.diag(np.ones(j), -r + j)
        else:
            Ap[:r, j - r : j] += inter[j] * np.diag(np.ones(r), 0)
    BB = np.zeros((r, m))
    BB[:r, : len(boundary_condition[0])] = boundary_condition

    def bn_func(t):
        return A.dot(gn(t))

    return A.dot(BB) + Ap, bn_func


def boundary_quasi_toep_to_boundary(quasi_toep_boundary, bn, inter, center):
    """Transformation from a boundary condition written as the extraction of the r first rows of the Quasi-Toeplitz matrix to a boundary condition written with ghost points.
    :param quasi_toep_boundary: Boundary condition BB written as the extraction of the r first rows of a Quasi-Toeplitz matrix
    :type quasi_toep_boundary: numpy.ndarray
    :param bn: Boundary data bn such that (U_0^{n+1},...,U_{r-1}^{n+1}) = BB (U_0^n,...,U_{m-1}^n) + b_n(t^n)
    :type bn: function
    :param inter: List of float coefficient that represent the interior scheme
    :type inter: list
    :param center: Index of the central coefficient of the scheme
    :type center: int
    :return: Boundary condition B written as U_{-r} = ..., ..., U_{-1} =...
    and the function g_n such that (U_{-r}^n,...,U_{-1}^n) = B(U_0^n,...,U_{m-1}^n) + g_n(t^n)
    :rtype: (numpy.ndarray,function)
    """
    r = center
    p = len(inter) - r - 1
    m = max(len(quasi_toep_boundary[0]), r + p)
    A = np.zeros((r, r))
    for j in range(r):
        A += inter[j] * np.diag(np.ones(r - j), j)
    Ap = np.zeros((r, m))
    for j in range(1, r + p + 1):
        if j < r:
            Ap[:r, :r] += inter[j] * np.diag(np.ones(j), -r + j)
        else:
            Ap[:r, j - r : j] += inter[j] * np.diag(np.ones(r), 0)
    BB = np.zeros((r, m))
    BB[:r, : len(quasi_toep_boundary[0])] = quasi_toep_boundary

    def gn_func(t):
        return np.linalg.inv(A).dot(bn(t))

    return np.linalg.inv(A).dot(BB - Ap), gn_func


def recflux(N, lamb):
    """Computes the coefficient of the Lax-Wendroff scheme at any order. Recursive function.

    :param N: Order of Lax Wendroff scheme
    :type N: int
    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :return: Parameters of the Lax-Wendroff scheme : interior coefficient and the center.
    :rtype: (np.ndarray,int)
    """
    if N == 1:
        return np.array([1]), 0
    else:
        m = N // 2
        M = (N - 1) // 2
        a = np.zeros(N)
        for k in range(N):
            a[N - 1 - k] = (-1) ** (k + N) * coefBinomial(N - 1, k)
        center = N - 1 - m
        prod = 1
        for i in range(-m, M + 1):
            if i != 0:
                prod *= lamb + i
        a = -a * prod / factorial(N)
        return schemsum(recflux(N - 1, lamb), (a, center))


def Lwconstructor(order):
    """Constructor of Lax-Wendroff scheme.

    :param order: Order of Lax Wendroff scheme
    :type order: int
    :return: Function depending on the Courant number lambda and maps to the parameters of the Lax-Wendroff scheme
    :rtype: function
    """

    def LWconstructorlamb(lamb):
        fluxschem = recflux(order, lamb)
        newfluxschem = -fluxschem[0], fluxschem[1] + 1
        susfluxschem = schemscal(schemsum(fluxschem, newfluxschem), -lamb)
        return schemsum((np.array([1]), 0), susfluxschem)

    return LWconstructorlamb
