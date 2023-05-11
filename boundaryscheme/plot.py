"""
This file defines utils function specific to the library
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt

from utils import boundary_to_boundary_quasi_toep, sort, epsilon, eta_func, parametrization
from boundaries import Dirichlet


def detKLplotsimple(s, n_param = 300, parametrization_bool = True):
    """
    draw the Kreiss-Lopatinskii curve for a numerical scheme s of the class Scheme 
    """
    fDKL = s.DKL()
    Param, Det = s.detKL(n_param, fDKL, parametrization_bool)
    plt.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2)
    plt.axvline(x=0, color="0.5")
    plt.axhline(y=0, color="0.5")
    plt.axis('equal')
    plt.show()


def detKLplot(schem, left_bound = Dirichlet(), lamb = None, sigma = 0, lambdacursor = False, sigmacursor = False):
    """
    draw the Kreiss-Lopatinskii curve for different lambdas and different sigmas or with cursor(s) 

    arguments : 
    schem :  a numerical scheme of the class Scheme depending on a parameter lambda
    left_bound : a left boundary of the class Boundary
    lamb : a value lambda or a numpy.ndarray of several values of lambda
    sigma : a real parameter representing the gap between the mesh and the boundary condition or a numpy.ndarray of several values of sigma
    lambdacursor : boolean to indicate the use of a cursor for lambda's moving among the lamb values
    sigmacursor : boolean to indicate the use of a cursor for sigma's moving among the sigma values
    """
    if lamb == None:
        lamb = np.array([schem(1).CFL])
    if isinstance(lamb, (int,float)):
        lamb = np.array([lamb], dtype=float)
    if isinstance(sigma, (int,float)):
        sigma = np.array([sigma], dtype=float)
    if sigmacursor == False:
        if lambdacursor == False:
            if len(sigma)==1:
                for l in lamb:
                    pass
            elif len(lamb) == 1:
                for sig in sigma:
                    pass
            else:
                for l in lamb:
                    for sig in sigma:
                        pass
        else:
            lambmax = np.max(lamb)
            if len(lamb) == 1 or np.min(lamb) < 0.001:
                lambmin = 0.001
            else:
                lambmin = np.min(lamb)
            pass
    else:
        if len(sigma) == 1:
            sigmax = 1/2
            sigmin = -1/2
        else:
            sigmax = np.max(sigma)
            sigmin = np.min(sigma)
        if lambdacursor == False:
            pass

        else:
            lambmax = np.max(lamb)
            if len(lamb) == 1 or np.min(lamb) < 0.001:
                lambmin = 0.001
            else:
                lambmin = np.min(lamb)
            pass
        






def symbolplot():
    pass

def numberzerosdetKL():
    pass