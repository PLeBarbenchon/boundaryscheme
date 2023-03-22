"""
This file defines the schemes classes
"""

import copy
import sympy as sp
import numpy as np
from math import *
from cmath import *
import numpy.polynomial.polynomial as nppol

from utils import boundary_to_boundary_quasi_toep, sort, epsilon, eta_func, parametrization
from boundaries import Dirichlet


class Scheme:
    def __init__(self, int, center, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        t = kwargs.get("t", 0)
        J = kwargs.get("J", 50)
        dx = kwargs.get("dx", 0.1)
        a = kwargs.get("a", 1)
        g = kwargs.get("g", lambda x: 0)
        dg = kwargs.get("dg", lambda x: 0)
        d2g = kwargs.get("d2g", lambda x: 0)
        d3g = kwargs.get("d3g", lambda x: 0)
        self.int = int
        self.center = center
        self.r = center
        self.p = len(self.int) - self.r - 1
        self.boundaryname = boundary.name()
        littleboundary, self.gn = boundary(center, sigma=self.sigma, dx=dx, a=a, g=g, dg=dg, d2g=d2g, d3g=d3g)
        self.m = max(len(littleboundary[0]), self.r + self.p)
        self.boundary = np.zeros((self.r, self.m))
        self.boundary[: self.r, : len(littleboundary[0])] = littleboundary
        self.boundary_quasi_toep, self.bn = boundary_to_boundary_quasi_toep(self.boundary, self.gn, int, center)

    def scheme(self):
        return self.int, self.center

    def toep(self, J, Borddroit=np.array([[]])):
        A = np.zeros((J, J))
        for k in range(len(self.int)):
            A += self.int[k] * np.diag(np.ones(J - abs(k - self.center)), k - self.center)
        A[0 : len(self.boundary_quasi_toep), 0 : len(self.boundary_quasi_toep[0])] = self.boundary_quasi_toep
        A[J - len(Borddroit) : J, J - len(Borddroit[0]) : J] = Borddroit
        return A

    def b_n(self, t, J):
        bn_vect = np.zeros(J)
        bn_vect[: self.center] = self.bn(t)
        return bn_vect

    def toep_circ(self, J):
        C = np.zeros((J, J))
        for k in range(len(self.int)):
            C += self.int[k] * np.diag(np.ones(J - abs(k - self.center)), k - self.center)
            if k - self.center >= 0:
                C += self.int[k] * np.diag(np.ones(abs(k - self.center)), -J + k - self.center)
            else:
                C += self.int[k] * np.diag(np.ones(abs(k - self.center)), J + k - self.center)
        return np.matrix(C)

    def symbol(self, n=300):
        T = np.linspace(0, 7, n)
        X = np.zeros_like(T)
        Y = np.zeros_like(T)
        for i in range(len(self.int)):
            X += self.int[i] * np.cos((i - self.center) * T)
            Y += self.int[i] * np.sin((i - self.center) * T)
        return X, Y

    def pol(self, z):
        r = self.center
        monome = nppol.Polynomial([0, 1])
        P = nppol.Polynomial(self.int) - z * monome**r
        return P

    def roots(self, z):
        Racinestotales = self.pol(z).roots()
        return sort(Racinestotales)

    def count_root(self, eta, eps, z0, kappa):
        z = z0 + eta * z0 / (2 * abs(z0))
        NewRoots = self.roots(z)
        selection = list(filter(lambda k: abs(k - kappa) < eps, NewRoots))
        return len(list(filter(lambda k: abs(k) < 1, selection)))

    def Kappa(self, z0):
        """
        selection of kappas that come form the inside of the unit disk
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

    def DKL(self): 
        z0 = sp.Symbol("z0", imaginary = True)
        B = copy.deepcopy(self.boundary_quasi_toep)
        B = sp.Matrix(B) - z0*sp.eye(self.r,self.m)
        b = sp.ones(1,self.r+1)
        for i in range (self.r+1):
            b[i] = sp.Symbol("b"+str(i), imaginary = True)
        for j in range (self.m-self.r):
            row = sp.zeros(1,self.m)
            for k in range (self.r+1):
                row[self.m-self.r-j+k-1] = b[k]
            B = B - np.dot(B[:,self.m-1-j],row)
        fdetKL = sp.lambdify([z0,b],sp.det(B[:self.r,:self.r]),"numpy")
        return fdetKL

    def detKL(self, n_param, DKL_formula, parametrization_bool):
        def scalar_detKL(self,z):
            Rz = nppol.polyfromroots(self.Kappa(z))
            return DKL_formula(z,Rz)
        if parametrization_bool:
            return parametrization(n_param,lambda z:scalar_detKL(self,z))
        else:
            param = np.linspace(0,2*pi,n_param)
            return param, np.vectorize(scalar_detKL)(self,np.exp(1j*param))

    def shortname(self):
        r = self.center
        s = "U_j^{n+1} = " + str(round(self.int[0], 3)) + " U_{j-" + str(r) + "}^n + "
        for i in range(1, len(self.int) - 1):
            if r > i:
                s += str(round(self.int[i], 3)) + " U_{j-" + str(r - i) + "}^n + "
            elif r == i:
                s += str(round(self.int[i], 3)) + " U_j^n + "
            else:
                s += str(round(self.int[i], 3)) + " U_{j+" + str(i - r) + "}^n + "
        if r > len(self.int) - 1:
            s += str(round(self.int[-1], 3)) + " U_{j-" + str(r - i) + "}^n"
        elif r == len(self.int) - 1:
            s += str(round(self.int[-1], 3)) + " U_j^n"
        else:
            s += str(round(self.int[-1], 3)) + " U_{j+" + str(i - r + 1) + "}^n"
        return s

    def name(self, boundary_bool=False, sigma_bool=False, lambda_bool=False):
        string = self.shortname()
        if boundary_bool:
            string += r" with boundary " + self.boundaryname
        if sigma_bool and self.sigma != 0:
            string += r" with $\sigma$ = " +str(self.sigma)
        if lambda_bool:
            string += r" for $\lambda$ = " + str(round(self.lamb,2))
        return string
    

    def __repr__(self):
        r = self.center
        s = "the interior equation of the scheme is :\n" + self.shortname() + "\nthe boundary equation is :"
        for i in range(len(self.boundary)):
            s += "\nU_{-" + str(r - i) + "}^n = "
            for j in range(len(self.boundary[i]) - 1):
                s += str(round(self.boundary[i, j], 3)) + " U_" + str(j) + "^n + "
            s += str(round(self.boundary[i, -1], 3)) + " U_" + str(len(self.boundary[0]) - 1) + "^n"
        return s


# rajouter une méthode pour plotter les solutions en fonction de x et t
# rajouter une méthode pour plotter les solutions en fonction de x avec animation en temps


class SchemeP0(Scheme):
    def __init__(self, int, center, boundary, **kwargs):
        super().__init__(int=int, center=center, boundary=boundary, **kwargs)
        assert len(int) == center + 1

    def detKL(self,n_param,DKL_formula, parametrization_bool):
        def scalar_detKL(self,z):
            b = np.array([self.int[i]/(self.int[-1]-z) for i in range (len(self.int))])
            b[-1] = 1
            return DKL_formula(z,b)
        if parametrization_bool:
            return parametrization(n_param,lambda z:scalar_detKL(self,z))
        else:
            param = np.linspace(0,2*pi,n_param)
            return param, np.vectorize(scalar_detKL)(self,np.exp(1j*param))



class O3Upwind4(SchemeP0):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        a0 = kwargs.get("a0", 0)
        a5 = kwargs.get("a5", 0)
        a6 = kwargs.get("a6", 0)
        self.lamb = lamb
        self.a0 = a0
        self.a5 = a5
        self.a6 = a6
        a1 = -(lamb**3) / 6 + 3 * lamb**2 / 2 - 13 * lamb / 3 + a5 + 4 * a6 - 4 * a0 + 4
        a2 = (lamb**3 - 8 * lamb**2 + 19 * lamb - 2 * (4 * a5 + 15 * a6 - 6 * a0 + 6)) / 2
        a3 = -(lamb**3) / 2 + 7 * lamb**2 / 2 - 7 * lamb + 6 * a5 + 20 * a6 - 4 * a0 + 4
        a4 = lamb**3 / 6 - lamb**2 + 11 * lamb / 6 - 4 * a5 - 10 * a6 + a0 - 1
        self.int = np.array([a6, a5, a4, a3, a2, a1, a0])
        self.center = 6
        while self.int[0] == 0:
            self.int = self.int[1:]
            self.center -= 1
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def shortname(self):
        return "O3Upwind"



class BWUpwind(SchemeP0):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        theta = kwargs.get("theta", 1 / 2)
        self.theta = theta
        self.lamb = lamb
        coeff1 = theta * (lamb - 1) * (lamb - 2) / 2 + (1 - theta) * (1 - lamb)
        coeff2 = theta * lamb * (2 - lamb) + (1 - theta) * lamb
        coeff3 = theta * lamb * (lamb - 1) / 2
        self.int = np.array([coeff3, coeff2, coeff1])
        self.center = 2
        while self.int[0] == 0:
            self.int = self.int[1:]
            self.center -= 1
        super().__init__(**kwargs)


class BeamWarming(SchemeP0):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        self.p = 0
        self.r = 2
        coeff1 = (lamb - 1) * (lamb - 2) / 2
        coeff2 = lamb * (2 - lamb)
        coeff3 = lamb * (lamb - 1) / 2
        self.int = np.array([coeff3, coeff2, coeff1])
        self.center = 2
        while self.int[0] == 0:
            self.int = self.int[1:]
            self.center -= 1
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 2

    def shortname(self):
        return "BeamWarming"



class Upwind(SchemeP0):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        self.p = 0
        self.r = 1
        coeff1 = 1 - lamb
        coeff2 = lamb
        self.int = np.array([coeff2, coeff1])
        self.center = 1
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return "Upwind"



class LaxWendroff(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        coeff1 = -(lamb - lamb**2) / 2
        coeff2 = 1 - lamb**2
        coeff3 = (lamb**2 + lamb) / 2
        self.int = np.array([coeff3, coeff2, coeff1])
        self.center = 1
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return "LaxWendroff"



class LaxFriedrichs(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        self.int = np.array([(1 + lamb) / 2, 0, (1 - lamb) / 2])
        self.center = 1
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return "LaxFriedrichs"



class ThirdOrder(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        self.r = 2
        self.p = 1
        self.int = np.array(
            [
                lamb**3 / 6 - lamb / 6,
                lamb + lamb * lamb / 2 - lamb**3 / 2,
                1 - lamb / 2 - lamb * lamb + lamb**3 / 2,
                -(lamb**3) / 6 + lamb * lamb / 2 - lamb / 3,
            ]
        )
        self.center = 2
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return "ThirdOrder"



class BB(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.lamb = lamb
        self.int = np.array([self.lamb / 4, self.lamb / 4, 1 - self.lamb / 4, -self.lamb / 4])
        self.center = 2
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return "BB"



class Dissipatif(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.D = kwargs.get("D", 0)
        self.lamb = lamb
        l = 2
        self.int = np.array([lamb / (2 * l), self.D, 1 - 2 * self.D, self.D, -lamb / (2 * l)])
        self.center = 2
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def shortname(self):
        return "Dissipatif"



class LW(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), **kwargs):
        self.sigma = kwargs.get("sigma", 0)
        self.order = kwargs.get("order", 2)
        self.lamb = lamb
        if self.order == 2:
            coeff1 = -(lamb - lamb**2) / 2
            coeff2 = 1 - lamb**2
            coeff3 = (lamb**2 + lamb) / 2
            self.int = np.array([coeff3, coeff2, coeff1])
            self.center = 1
        elif self.order == 3:
            self.int = np.array(
                [
                    -lamb * (1 - lamb**2) / 6,
                    -lamb * (lamb + 1) * (lamb - 2) / 2,
                    1 - lamb * (1 + 2 * lamb - lamb**2) / 2,
                    -lamb * (lamb - 1) * (lamb - 2) / 6,
                ]
            )
            self.center = 2
        elif self.order == 4:
            self.int = np.array(
                [
                    lamb * (lamb - 1) * (lamb + 1) * (lamb + 2) / 24,
                    -lamb * (lamb - 2) * (lamb + 2) * (lamb + 1) / 6,
                    1 + lamb**2 * (lamb**2 - 5) / 4,
                    -lamb * (lamb - 2) * (lamb + 2) * (lamb - 1) / 6,
                    lamb * (lamb - 1) * (lamb + 1) * (lamb - 2) / 24,
                ]
            )
            self.center = 2
        elif self.order == 5:
            self.int = np.array(
                [
                    lamb * (lamb - 2) * (lamb - 1) * (lamb + 1) * (lamb + 2) / 120,
                    -lamb * (lamb - 1) * (lamb - 3) * (lamb + 1) * (lamb + 2) / 24,
                    lamb * (lamb - 2) * (lamb - 3) * (lamb + 1) * (lamb + 2) / 12,
                    1 - lamb * (lamb**4 - 3 * lamb**3 - 5 * lamb**2 + 15 * lamb + 4) / 12,
                    lamb * (lamb - 1) * (lamb - 2) * (lamb - 3) * (lamb + 2) / 24,
                    -lamb * (lamb - 1) * (lamb - 2) * (lamb - 3) * (lamb + 1) / 120,
                ]
            )
            self.center = 3
        elif self.order == 6:
            self.int = np.array(
                [
                    lamb * (lamb - 2) * (lamb + 3) * (lamb - 1) * (lamb + 1) * (lamb + 2) / 720,
                    -lamb * (lamb - 1) * (lamb - 3) * (lamb + 1) * (lamb + 2) * (lamb + 3) / 120,
                    lamb * (lamb - 2) * (lamb - 3) * (lamb + 1) * (lamb + 3) * (lamb + 2) / 48,
                    1 - lamb * lamb * (lamb**2 - 7) * (lamb**2 - 7) / 36,
                    lamb * (lamb - 1) * (lamb - 2) * (lamb - 3) * (lamb + 3) * (lamb + 2) / 48,
                    -lamb * (lamb - 1) * (lamb - 2) * (lamb - 3) * (lamb + 1) * (lamb + 3) / 120,
                    lamb * (lamb - 1) * (lamb - 2) * (lamb - 3) * (lamb + 2) * (lamb + 1) / 720,
                ]
            )
            self.center = 3
        super().__init__(int=self.int, center=self.center, boundary=boundary, **kwargs)

    def CFL(self):
        return 1

    def shortname(self):
        return f"Lax Wendroff {self.order}"



class LeapFrog(Scheme):
    pass


class CrankNicholson(Scheme):
    pass
