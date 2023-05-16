"""
This file defines utils function specific to the library
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from .utils import boundary_to_boundary_quasi_toep, sort, epsilon, eta_func, parametrization
from .boundaries import Dirichlet, SILW
from .schemes import *


def detKLplotsimple(s, n_param = 300, parametrization_bool = True):
    """
    draw the Kreiss-Lopatinskii curve for a numerical scheme s of the class Scheme 
    """
    Param, Det = s.detKL(n_param, parametrization_bool)
    plt.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2)
    plt.axvline(x=0, color="0.5")
    plt.axhline(y=0, color="0.5")
    plt.axis('equal')
    plt.show()

# detKLplotsimple(BeamWarming(1, SILW(2,3)),parametrization_bool=False)


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
        






def symbolplot(schem, lamb = None, order = 2, lambdacursor = False, nparam=300):
    if lambdacursor == False:
        fig = plt.figure(1, figsize=(10, 8))
        if isinstance(lamb, (int,float)):
            S = schem(lamb,order = order)
            X,Y = S.symbol(nparam)
            T = np.linspace(0,7,1000)

            plt.ylim(-1,1)
            plt.xlim(-1,1)
            plt.plot(np.cos(T),np.sin(T),"--",  color = "black", label = "unit circle")
            plt.plot(X,Y, linewidth=2, label = "symbol")
            plt.axis("equal")
            plt.legend(loc="best")
            plt.title("Symbol of "+ S.name(lambda_bool = True))
        else:
            assert (type(lamb) == np.ndarray and lamb.ndim == 1) or type(lamb) == list
            for l in lamb:
                S = schem(l,order = order)
                X,Y = S.symbol(nparam)

                plt.ylim(-1,1)
                plt.xlim(-1,1)
                plt.plot(X,Y, linewidth=2, label = f"$\lambda = $ {l}")
            T = np.linspace(0,7,1000)
            plt.plot(np.cos(T),np.sin(T),"--", color = "black", label = "unit circle")
            plt.axis("equal")
            plt.legend(loc="upper left")
            plt.title("Symbol of "+ S.name(lambda_bool = False))

    else:
        if lamb == None:
            lamb = schem(1).CFL
        if isinstance(lamb, (int,float)):
            lamb = np.array([lamb], dtype=float)
        fig = plt.figure(1, figsize=(10, 8))
        plt.title("Symbol of "+ schem(1).name(lambda_bool = False))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.25, bottom=0.25)
        ax.set_xlim(-3,3)
        ax.axis("equal")

        if len(lamb) == 1:
            lamb0 = lamb[0]/2
        else:
            lamb0 = (np.min(lamb) +np.max(lamb))/2

        Param = np.linspace(0,1,nparam)
        Param = np.cos(2*pi*Param) + 1j * np.sin(2*pi*Param)
        Theta = np.linspace(0,2*pi,500)

        X,Y = schem(lamb0, order = order).symbol(nparam)
        [line] = ax.plot(X,Y, linewidth=2, label = "symbol")
        [cible] = ax.plot(np.cos(Theta),np.sin(Theta),"--",  color = "black", label = "unit circle")
        ax.legend(loc="best")
        ax.set_xlim(-6, 2)
        ax.axis("equal")
        if len(lamb) == 1:
            minvalue = 0.0001
            maxvalue = lamb[0]
        else:
            minvalue = np.min(lamb)
            maxvalue = np.max(lamb)

        axlamb = plt.axes([0.25, 0.08, 0.65, 0.03])
        lambda_slider = Slider(
            ax=axlamb,
            label="lambda",
            valmin=minvalue,
            valmax=maxvalue,
            valinit=lamb0,
        )

    

        def update(val):
            S = schem(lambda_slider.val,order=order)
            X,Y = S.symbol(nparam)
            line.set_xdata(X)
            line.set_ydata(Y)
            fig.canvas.draw_idle()

        # register the update function with each slider
        lambda_slider.on_changed(update)

        # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', hovercolor='0.975')

        def reset(event):
            lambda_slider.reset()
        button.on_clicked(reset)
    plt.show()





def numberzerosdetKL():
    pass

