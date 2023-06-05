"""
This file defines the schemes classes
"""

import copy
import sympy as sp
import numpy as np
from math import *
from cmath import *
import numpy.polynomial.polynomial as nppol

from .utils import boundary_to_boundary_quasi_toep, sort, epsilon, eta_func, parametrization, Lwconstructor
from .boundaries import Dirichlet


class Scheme:
    """This is a class to represent a numerical scheme.

    :param inter: List of float coefficient that represent the interior scheme
    :type inter: list
    :param center: Index of the central coefficient of the scheme
    :type center: int
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional

    .. seealso:: [Boutin, Le Barbenchon, Seguin, 2023 : Stability of finite difference schemes for the hyperbolic initial boundary value problem by winding number computations]
    """

    def __init__(self, inter, center, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.inter = inter
        self.center = center
        self.r = center
        self.p = len(self.inter) - self.r - 1
        self.boundaryname = boundary.name()
        littleboundary, self.gn = boundary(center, sigma=self.sigma)
        self.m = max(len(littleboundary[0]), self.r + self.p)
        self.boundary = np.zeros((self.r, self.m))
        self.boundary[: self.r, : len(littleboundary[0])] = littleboundary
        self.boundary_quasi_toep, self.bn = boundary_to_boundary_quasi_toep(self.boundary, self.gn, inter, center)

        z0 = sp.Symbol("z0", imaginary=True)
        B = copy.deepcopy(self.boundary_quasi_toep)
        B = sp.Matrix(B) - z0 * sp.eye(self.r, self.m)
        b = sp.ones(1, self.r + 1)
        for i in range(self.r + 1):
            b[i] = sp.Symbol("b" + str(i), imaginary=True)
        for j in range(self.m - self.r):
            row = sp.zeros(1, self.m)
            for k in range(self.r + 1):
                row[self.m - self.r - j + k - 1] = b[k]
            B = B - np.dot(B[:, self.m - 1 - j], row)
        self.DKL = sp.lambdify([z0, b], sp.det(B[: self.r, : self.r]), "numpy")

    def scheme(self):
        """Returns the parameters of the scheme : the interior coefficients et the center
        :return: Parameters of the scheme
        :rtype: (list,int)
        """
        return self.inter, self.center

    def toep(self, J, right_bound=np.array([[]])):
        """Returns the Toeplitz matrix of the scheme with two boundaries
        :param J: Size of the Toeplitz matrix
        :type J: int
        :param right_bound: Right boundary of the scheme
        :type right_bound: class:`Boundary`
        :return: Toeplitz matrix of the scheme
        :rtype: numpy.ndarray
        """
        A = np.zeros((J, J))
        for k in range(len(self.inter)):
            A += self.inter[k] * np.diag(np.ones(J - abs(k - self.center)), k - self.center)
        A[0 : len(self.boundary_quasi_toep), 0 : len(self.boundary_quasi_toep[0])] = self.boundary_quasi_toep
        A[J - len(right_bound) : J, J - len(right_bound[0]) : J] = right_bound
        return A

    def toep_circ(self, J):
        """Returns the circulant Toeplitz matrix of the scheme
        :param J: Size of the circulant matrix
        :type J: int
        :return: Circulant matrix of the scheme
        :rtype: numpy.ndarray
        """
        C = np.zeros((J, J))
        for k in range(len(self.inter)):
            C += self.inter[k] * np.diag(np.ones(J - abs(k - self.center)), k - self.center)
            if k - self.center >= 0:
                C += self.inter[k] * np.diag(np.ones(abs(k - self.center)), -J + k - self.center)
            else:
                C += self.inter[k] * np.diag(np.ones(abs(k - self.center)), J + k - self.center)
        return np.matrix(C)

    def symbol(self, n=300):
        """To draw the symbol curve of the scheme
        :param n: Number of points of the parametrisation of the unit circle, defaults to 300
        :type n: int, optional
        :return: Abscissa and Ordinate of the symbol curve of the scheme
        :rtype: (list,list)
        """
        T = np.linspace(0, 7, n)
        X = np.zeros_like(T)
        Y = np.zeros_like(T)
        for i in range(len(self.inter)):
            X += self.inter[i] * np.cos((i - self.center) * T)
            Y += self.inter[i] * np.sin((i - self.center) * T)
        return X, Y

    def pol(self, z):
        """Returns the characteristic polynomial of the scheme
        :param z: Parameter z of the characteristic polynomial
        :type z: complex
        :return: Characteristic polynomial of the scheme
        :rtype: class: Polynomial
        """
        r = self.center
        monome = nppol.Polynomial([0, 1])
        P = nppol.Polynomial(self.inter) - z * monome**r
        return P

    def roots(self, z):
        """Returns the roots of the characteristic polynomial
        :param z: Parameter z of the characteristic polynomial
        :type z: complex
        :return: Roots of the characteristic polynomial sorting by utils.sort function
        :rtype: list
        """
        Racinestotales = self.pol(z).roots()
        return sort(Racinestotales)

    def count_root(self, eta, eps, z0, kappa):
        """Counts the number of roots from kappa(z0) that comes from the inside of the unit disk
        :param eta: Gap between z and z0
        :type eta: float
        :param eps: Threshold
        :type eps: float
        :param z0: Parameter z0
        :type z0: complex
        :param kappa: Multiple root kappa linked to z0
        :type kappa: complex
        :return: Number of roots coming from the inside of the unit disk
        :rtype: int
        """
        z = z0 + eta * z0 / (2 * abs(z0))
        NewRoots = self.roots(z)
        selection = list(filter(lambda k: abs(k - kappa) < eps, NewRoots))
        return len(list(filter(lambda k: abs(k) < 1, selection)))

    def Kappa(self, z0):
        """Selection of kappas that come from the inside of the unit disk
        :param z0: Parameter z0 of the characteristic polynomial
        :type z0: complex
        :return: Roots of the characteristic polynomial for z0 that are inside the unit disk or are coming from the inside of the unit disk
        :rtype: list
        """
        delta = 10 ** (-10)
        Racinestotales = self.roots(z0)
        eps = epsilon(Racinestotales)
        RootsFromInside = []
        stock = []
        for x in Racinestotales:
            if abs(x) < 1 - delta:
                RootsFromInside.append(x)
            elif abs(x) < 1 + delta and x not in stock:
                stock.append(x)
                eta = eta_func(eps, x, 1000, self.pol(z0), self.center)
                n = self.count_root(eta, eps, z0, x)
                for i in range(n):
                    RootsFromInside.append(x)
        assert len(RootsFromInside) == self.center
        return RootsFromInside

    def detKL(self, n_param, parametrization_bool):
        """Computes the Kreiss-Lopatinskii determinant on the unit circle
        :param n_param: Number of discretization of the unit circle
        :type n_param: int
        :param parametrization_bool: If True, the discretization is refined close to the origin. If False, the discretization is uniform on the unit circle.
        :type parametrization_bool: bool
        :return: The parametrization of the unit circle and the value of the Kreiss-Lopatinskii determinant on the unit circle
        :rtype: (numpy.ndarray, numpy.ndarray)
        """

        def scalar_detKL(self, z):
            Rz = nppol.polyfromroots(self.Kappa(z))
            return self.DKL(z, Rz)

        if parametrization_bool:
            return parametrization(n_param, lambda z: scalar_detKL(self, z))
        else:
            param = np.linspace(0, 2 * pi, n_param)
            return param, np.vectorize(scalar_detKL)(self, np.exp(1j * param))

    def shortname(self):
        """Gives a short name of the scheme"""
        r = self.center
        s = "U_j^{n+1} = " + str(round(self.inter[0], 3)) + " U_{j-" + str(r) + "}^n + "
        for i in range(1, len(self.inter) - 1):
            if r > i:
                s += str(round(self.inter[i], 3)) + " U_{j-" + str(r - i) + "}^n + "
            elif r == i:
                s += str(round(self.inter[i], 3)) + " U_j^n + "
            else:
                s += str(round(self.inter[i], 3)) + " U_{j+" + str(i - r) + "}^n + "
        if r > len(self.inter) - 1:
            s += str(round(self.inter[-1], 3)) + " U_{j-" + str(r - i) + "}^n"
        elif r == len(self.inter) - 1:
            s += str(round(self.inter[-1], 3)) + " U_j^n"
        else:
            s += str(round(self.inter[-1], 3)) + " U_{j+" + str(i - r + 1) + "}^n"
        return s

    def name(self, boundary_bool=False, sigma_bool=False, lambda_bool=False):
        """Gives a name of the scheme"""
        string = self.shortname()
        if boundary_bool:
            string += r" with boundary " + self.boundaryname
        if sigma_bool and self.sigma != 0:
            string += r" with $\sigma$ = " + str(self.sigma)
        if lambda_bool:
            string += r" for $\lambda$ = " + str(round(self.lamb, 2))
        return string

    def __repr__(self):
        """Print the numerical scheme"""
        r = self.center
        s = "the interior equation of the scheme is :\n" + self.shortname() + "\nthe boundary equation is :"
        for i in range(len(self.boundary)):
            s += "\nU_{-" + str(r - i) + "}^n = "
            for j in range(len(self.boundary[i]) - 1):
                s += str(round(self.boundary[i, j], 3)) + " U_" + str(j) + "^n + "
            s += str(round(self.boundary[i, -1], 3)) + " U_" + str(len(self.boundary[0]) - 1) + "^n"
        return s


class SchemeP0(Scheme):
    """This is a class to represent a totally upwind scheme (with p = 0).

    :param inter: List of float coefficient that represent the interior scheme
    :type inter: list
    :param center: Index of the central coefficient of the scheme
    :type center: int
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional

    .. seealso:: [Boutin, Le Barbenchon, Seguin, 2023 : On the stability of totally upwind schemes for the hyperbolic initial boundary value problem]
    """

    def __init__(self, inter, center, boundary, sigma=0):
        """Constructor Method"""
        super().__init__(inter=inter, center=center, boundary=boundary, sigma=sigma)
        assert len(inter) == center + 1

    def detKL(self, n_param, parametrization_bool):
        """
        Redefinition of the detKL method of the :class: `Scheme` but in the particular case of totally upwind scheme
        :param n_param: Number of discretization of the unit circle
        :type n_param: int
        :param parametrization_bool: If True, the discretization is refined close to the origin. If False, the discretization is uniform on the unit circle.
        :type parametrization_bool: bool
        :return: The parametrization of the unit circle and the value of the Kreiss-Lopatinskii determinant on the unit circle
        :rtype: (numpy.ndarray, numpy.ndarray)
        """

        def scalar_detKL(self, z):
            b = np.array([self.inter[i] / (self.inter[-1] - z) for i in range(len(self.inter))])
            b[-1] = 1
            return self.DKL(z, b)

        if parametrization_bool:
            return parametrization(n_param, lambda z: scalar_detKL(self, z))
        else:
            param = np.linspace(0, 2 * pi, n_param)
            return param, np.vectorize(scalar_detKL)(self, np.exp(1j * param))


class BeamWarming(SchemeP0):
    """This is a class to represent the Beam-Warming scheme (which is totally upwind).

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.p = 0
        self.r = 2
        coeff1 = (lamb - 1) * (lamb - 2) / 2
        coeff2 = lamb * (2 - lamb)
        coeff3 = lamb * (lamb - 1) / 2
        self.inter = np.array([coeff3, coeff2, coeff1])
        self.center = 2
        self.CFL = 2
        while self.inter[0] == 0:
            self.inter = self.inter[1:]
            self.center -= 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma)

    def shortname(self):
        """Name method"""
        return "BeamWarming"


class Upwind(SchemeP0):
    """This is a class to represent the Upwind scheme (which is totally upwind).

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.p = 0
        self.r = 1
        coeff1 = 1 - lamb
        coeff2 = lamb
        self.inter = np.array([coeff2, coeff1])
        self.center = 1
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "Upwind"


class LaxWendroff(Scheme):
    """This is a class to represent the Lax-Wendroff scheme (of order 2).

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        coeff1 = -(lamb - lamb**2) / 2
        coeff2 = 1 - lamb**2
        coeff3 = (lamb**2 + lamb) / 2
        self.inter = np.array([coeff3, coeff2, coeff1])
        self.center = 1
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "LaxWendroff"


class LaxFriedrichs(Scheme):
    """This is a class to represent the Lax-Friedrichs scheme.

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.inter = np.array([(1 + lamb) / 2, 0, (1 - lamb) / 2])
        self.center = 1
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "LaxFriedrichs"


class ThirdOrder(Scheme):
    """This is a class to represent the O3-scheme.

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.r = 2
        self.p = 1
        self.inter = np.array(
            [
                lamb**3 / 6 - lamb / 6,
                lamb + lamb * lamb / 2 - lamb**3 / 2,
                1 - lamb / 2 - lamb * lamb + lamb**3 / 2,
                -(lamb**3) / 6 + lamb * lamb / 2 - lamb / 3,
            ]
        )
        self.center = 2
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "ThirdOrder"


class BB(Scheme):
    """This is a class to represent an example of numerical scheme.

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.inter = np.array([self.lamb / 4, self.lamb / 4, 1 - self.lamb / 4, -self.lamb / 4])
        self.center = 2
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "BB"


class Dissipatif(Scheme):
    """This is a class to represent a dissipative scheme.

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    :param D: A parameter of the dissipativity, defaults to 0
    :type D: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, D=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.D = D
        self.lamb = lamb
        l = 2
        self.inter = np.array([lamb / (2 * l), self.D, 1 - 2 * self.D, self.D, -lamb / (2 * l)])
        self.center = 2
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "Dissipative"


class LW(Scheme):
    """This is a class to represent a Lax-Wendroff scheme of any order.

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    :param order: Order of the Lax-Wendroff scheme, defaults to 2
    :type order: int, optional
    """

    allordersschem = [None, None, Lwconstructor(2)]

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, order=2, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.order = order
        self.lamb = lamb
        n = len(LW.allordersschem) - 1
        if order > n:
            for i in range(order - n):
                LW.allordersschem.append(None)
            LW.allordersschem[order] = Lwconstructor(order)
        if LW.allordersschem[order] == None:
            LW.allordersschem[order] = Lwconstructor(order)
        self.inter, self.center = LW.allordersschem[order](lamb)
        self.CFL = 1
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return f"Lax Wendroff {self.order}"


# class LeapFrog(Scheme):
#     raise NotImplementedError('')


# class CrankNicholson(Scheme):
#     raise NotImplementedError('')
