"""
This file defines utils function specific to the library
"""

import numpy as np
from math import *
from cmath import *


def sort(L):
    """
    bubble sort
    sorting in the ascending order of modulus and for values with same modulus, sort in the ascending order of argument (between -pi and pi)
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
    """
    return the sector of z, it is an integer "a" between 0 and 7 which correspond to the angular sector [a*pi/4, (a+1)*pi/4[
    Warning : the argument is between -pi and pi and here "a" is between 0 and 7
    """
    ph = phase(z) * 4 / pi
    if ph < 0:
        return int(ph) + 7
    else:
        return int(ph)


def neighboor(sector1, z2):
    """
    sector1 is an integer between 0 and 7 and represent the angular sector [sector1*pi/4, (sector1+1)*pi/4[
    z2 is a complex number
    return True iff the complex z2 is in a neighboring sector of [sector1*pi/4, (sector1+1)*pi/4[
    """
    return abs((sector(z2) - sector1) % 8) <= 1



def parametrization(n_param, curve_formula):
    """
    n_param is an integer for the default discretization of [0,2pi]
    curve_formula is a fonction : z in the unit circle mapsto a complex and represent a curve
    return the discretization of the curve refine if it is needed (as [ZapataMartin2014] procedure).
    """
    dx = 2 * pi / n_param
    current_dx = dx
    current_param = 0
    Param = [current_param]
    curve = [curve_formula(1)]
    current_sector = sector(curve_formula(np.exp(1j * Param[0])))#calcul le secteur du premier point
    c = 0
    s = 0
    while Param[-1] < 2 * pi:
        current_param = Param[-1] + current_dx
        current_curve_point = curve_formula(np.exp(1j * current_param))
        if neighboor(current_sector,current_curve_point):
            Param.append(current_param)
            curve.append(current_curve_point)
            current_sector = sector(current_curve_point)
            current_dx = dx
            c = 0
        elif c<40:
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
    if Param[-1]>2*pi:
        Param[-1] = 0
        curve[-1] = curve_formula(1)
    return np.array(Param), np.array(curve)


def epsilon(L):
    """
    L is a list of elements x_i
    return min |x_i - x_j| / 2 (when x_i != x_j)
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
    """
    return min |polynom(kappa)|/(1+eps)^r for kappa on the circle centered in kappa0 of radius eps
    N is the number of discretization of the circle
    """
    theta = np.linspace(0, 2 * pi, N)
    circle = np.cos(theta) + 1j * np.sin(theta)
    kappas = kappa0 + eps * circle
    val_pol = np.abs(polynom(kappas))
    mini = min(val_pol)
    assert mini > 0
    return mini / ((1 + eps) ** r)


def boundary_to_boundary_quasi_toep(boundary_condition, gn, Int, center):
    """
    take a boundary condition written as U_{-r} = ..., ...,U_{-1} =..., etc
    return T_J the extraction of the r first rows of the Quasi-Toeplitz matrix and the function b_n(t) such that (U_0^{n+1},...,U_{r-1}^{n+1}) = T_J (U_0^n,...,U_{m-1}^n) + b_n(t^n)
    """
    r = center
    p = len(Int) - r - 1
    m = max(len(boundary_condition[0]), r + p)
    A = np.zeros((r, r))
    for j in range(r):
        A += Int[j] * np.diag(np.ones(r - j), j)
    Ap = np.zeros((r, m))
    for j in range(1, r + p + 1):
        if j < r:
            Ap[:r, :r] += Int[j] * np.diag(np.ones(j), -r + j)
        else:
            Ap[:r, j - r : j] += Int[j] * np.diag(np.ones(r), 0)
    BB = np.zeros((r, m))
    BB[:r, : len(boundary_condition[0])] = boundary_condition

    def bn_func(t):
        return A.dot(gn(t))

    return A.dot(BB) + Ap, bn_func


def boundary_quasi_toep_to_boundary(boundary, bn, Int, center):
    """
    take T_J (boundary) the extraction of the r first rows of a Quasi-Toeplitz matrix and a function bn such that (U_0^{n+1},...,U_{r-1}^{n+1}) = T_J (U_0^n,...,U_{m-1}^n) + b_n(t^n)
    return B the boundary condition written as U_{-r} = ..., ..., U_{-1} =...
    and the function g_n such that (U_{-r}^n,...,U_{-1}^n) = B(U_0^n,...,U_{m-1}^n) + g_n(t^n)
    """
    r = center
    p = len(Int) - r - 1
    m = max(len(boundary[0]), r + p)
    A = np.zeros((r, r))
    for j in range(r):
        A += Int[j] * np.diag(np.ones(r - j), j)
    Ap = np.zeros((r, m))
    for j in range(1, r + p + 1):
        if j < r:
            Ap[:r, :r] += Int[j] * np.diag(np.ones(j), -r + j)
        else:
            Ap[:r, j - r : j] += Int[j] * np.diag(np.ones(r), 0)
    BB = np.zeros((r, m))
    BB[:r, : len(boundary[0])] = boundary

    def gn_func(t):
        return np.linalg.inv(A).dot(bn(t))

    return np.linalg.inv(A).dot(BB - Ap), gn_func
