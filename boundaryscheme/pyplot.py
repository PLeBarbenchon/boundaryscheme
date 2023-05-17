"""
This file defines utils function specific to the library
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from .utils import boundary_to_boundary_quasi_toep, sort, epsilon, eta_func, parametrization
from .boundaries import Dirichlet
from .complex_winding_number import *


def show():
    plt.show()

def detKLplotsimple(s, n_param = 300, parametrization_bool = True):
    """
    draw the Kreiss-Lopatinskii curve for a numerical scheme s of the class Scheme 
    """
    ax = plt.subplot()
    Param, Det = s.detKL(n_param, parametrization_bool)
    ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2)
    ax.axvline(x=0, color="0.5")
    ax.axhline(y=0, color="0.5")
    ax.axis('equal')
    return ax



def detKLplot(schem, left_bound = Dirichlet(), lamb = None, sigma = 0, lambdacursor = False, sigmacursor = False, nparam = 300, parametrization_bool = True,  order = 2, fig_size = (6,4)):
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
    if lamb is None:
        lamb = np.array([schem(1).CFL])
    if isinstance(lamb, (int,float)):
        lamb = np.array([lamb], dtype=float)
    if isinstance(sigma, (int,float)):
        sigma = np.array([sigma], dtype=float)
    if sigmacursor == False:
        if lambdacursor == False:
            if len(sigma)==1:
                lamb = np.sort(lamb)
                fig,ax = plt.subplots(figsize=fig_size)
                for l in lamb:
                    S = schem(l, left_bound, order=order, sigma=sigma[0])
                    Param, Det = S.detKL(nparam, parametrization_bool)
                    ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, label=f"$\lambda$ = " + str(l))
                ax.axvline(x=0, color="0.5")
                ax.axhline(y=0, color="0.5")
                ax.legend(loc='best')
                plt.title("Symbol of "+ S.name(boundary_bool = True, sigma_bool = True,lambda_bool = False))
            elif len(lamb) == 1:
                sigma = np.sort(sigma)
                fig,ax = plt.subplots(figsize=fig_size)
                for sig in sigma:
                    S = schem(lamb[0], left_bound, order=order, sigma=sig)
                    Param, Det = S.detKL(nparam, parametrization_bool)
                    ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, label=f"$\sigma$ = " + str(sig))
                ax.axvline(x=0, color="0.5")
                ax.axhline(y=0, color="0.5")
                ax.legend(loc='best')
                plt.title("Symbol of "+ S.name(boundary_bool = True, sigma_bool = False,lambda_bool = True))
            else:
                sigma = np.sort(sigma)
                fig,ax = plt.subplots(figsize=fig_size)
                for l in lamb:
                    for sig in sigma:
                        S = schem(l, left_bound, order=order, sigma=sig)
                        Param, Det = S.detKL(nparam, parametrization_bool)
                        ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, label=f"$\lambda$ = " + str(l) + f"et $\sigma$ = " + str(sig))
                ax.axvline(x=0, color="0.5")
                ax.axhline(y=0, color="0.5")
                ax.legend(loc='best')
                plt.title("Symbol of "+ S.name(boundary_bool = True, sigma_bool = False,lambda_bool = False))
        else:
            if len(sigma)>1:
                raise TypeError("with a lambda cursor, sigma has to be a single value")
            lambmax = np.max(lamb)
            if len(lamb) == 1 or np.min(lamb) < 0.001:
                lambmin = 0.001
            else:
                lambmin = np.min(lamb)
            
            lamb0 = (lambmax + lambmin) / 2
            fig = plt.figure(1, figsize=(10, 8))
            plt.title("DKL curve of "+ schem(1, left_bound,sigma=sigma[0]).name(boundary_bool = True, sigma_bool=True))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.25, bottom=0.25)
            ax.set_xlim(-3,3)
            ax.axis("equal")


            def calc_det(l):
                S = schem(l, left_bound, sigma=sigma[0], order=order)
                return S.detKL(nparam, parametrization_bool)[1]

            Det = calc_det(lamb0)
            [line] = ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, color='red')
            ax.axvline(x=0, color="0.5")
            ax.axhline(y=0, color="0.5")
            ax.axis("equal")
            
            lambda_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
            lambda_slider = Slider(ax = lambda_slider_ax, label='lambda', valmin=lambmin, valmax=lambmax, valinit=lamb0)


        

            def update(val):
                Det = calc_det(lambda_slider.val)
                line.set_xdata([z.real for z in Det])
                line.set_ydata([z.imag for z in Det])
                xmin = min([z.real for z in Det])
                xmax = max([z.real for z in Det])
                ymin = min([z.imag for z in Det])
                ymax = max([z.imag for z in Det])
                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
                fig.canvas.draw_idle()

            # register the update function with each slider
            lambda_slider.on_changed(update)

            # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
            button = Button(resetax, 'Reset', hovercolor='0.975')

            def reset(event):
                lambda_slider.reset()
            button.on_clicked(reset)

    else:
        if len(sigma) == 1:
            sigmax = 1/2
            sigmin = -1/2
        else:
            sigmax = np.max(sigma)
            sigmin = np.min(sigma)
        if lambdacursor == False:
            if len(lamb)>1:
                raise TypeError("With a sigma cursor, lambda has to be a single value")
            sigma0 = (sigmax + sigmin) / 2
            fig = plt.figure(1, figsize=(10, 8))
            plt.title("DKL curve of "+ schem(lamb[0], left_bound).name(boundary_bool = True, lambda_bool = True))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.25, bottom=0.25)
            ax.set_xlim(-3,3)
            ax.axis("equal")


            def calc_det(sig):
                S = schem(lamb[0], left_bound, sigma=sig, order=order)
                return S.detKL(nparam, parametrization_bool)[1]

            Det = calc_det(sigma0)
            [line] = ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, color='red')
            ax.axvline(x=0, color="0.5")
            ax.axhline(y=0, color="0.5")
            ax.axis("equal")
            
            sigma_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
            sigma_slider = Slider(ax = sigma_slider_ax, label='sigma', valmin=sigmin, valmax=sigmax, valinit=sigma0)
        
            def update(val):
                Det = calc_det(sigma_slider.val)
                line.set_xdata([z.real for z in Det])
                line.set_ydata([z.imag for z in Det])
                xmin = min([z.real for z in Det])
                xmax = max([z.real for z in Det])
                ymin = min([z.imag for z in Det])
                ymax = max([z.imag for z in Det])
                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
                fig.canvas.draw_idle()

            # register the update function with each slider
            sigma_slider.on_changed(update)

            # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
            button = Button(resetax, 'Reset', hovercolor='0.975')

            def reset(event):
                sigma_slider.reset()
            button.on_clicked(reset)

        else:
            lambmax = np.max(lamb)
            if len(lamb) == 1 or np.min(lamb) < 0.001:
                lambmin = 0.001
            else:
                lambmin = np.min(lamb)
            
            sigma0 = (sigmax + sigmin)/2
            lamb0 = (lambmin + lambmax)/2

            fig = plt.figure(1, figsize=(10, 8))
            plt.title("DKL curve of "+ schem(1, left_bound).name(boundary_bool = True))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.25, bottom=0.25)
            ax.set_xlim(-3,3)
            ax.axis("equal")

            def calc_det(l, sigm):
                S = schem(l, left_bound, sigma=sigm, order=order)
                return S.detKL(nparam, parametrization_bool)[1]
            
            Det = calc_det(lamb0, sigma0)
            [line] = ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, color='red')
            plt.axvline(x=0, color="0.5")
            plt.axhline(y=0, color="0.5")


            axlambda = plt.axes([0.25, 0.1, 0.65, 0.03])
            lambda_slider = Slider(
                ax=axlambda,
                label='lambda',
                valmin=lambmin,
                valmax=lambmax,
                valinit=lamb0,
            )

            # Make a vertically oriented slider to control the amplitude
            axsigm = plt.axes([0.1, 0.25, 0.0225, 0.63])
            sigma_slider = Slider(
                ax=axsigm, label="sigma", valmin=sigmin, valmax=sigmax, valinit=sigma0, orientation="vertical"
            )

            def update(val):
                Det = calc_det(lambda_slider.val, sigma_slider.val)
                line.set_xdata([z.real for z in Det])
                line.set_ydata([z.imag for z in Det])
                xmin = min([z.real for z in Det])
                xmax = max([z.real for z in Det])
                ymin = min([z.imag for z in Det])
                ymax = max([z.imag for z in Det])
                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
                fig.canvas.draw_idle()
            
            lambda_slider.on_changed(update)
            sigma_slider.on_changed(update)
            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
            button = Button(resetax, 'Reset', hovercolor='0.975')

            def reset(event):
                lambda_slider.reset()
                sigma_slider.reset()

    button.on_clicked(reset)

    return ax
        






def symbolplot(schem, lamb = None, order = 2, lambdacursor = False, nparam=300, fig_size = (10,8)):
    if lambdacursor == False:
        fig,ax = plt.subplots(figsize=fig_size)
        if isinstance(lamb, (int,float)):
            S = schem(lamb,order = order)
            X,Y = S.symbol(nparam)
            T = np.linspace(0,7,1000)

            plt.ylim(-1,1)
            plt.xlim(-1,1)
            ax.plot(np.cos(T),np.sin(T),"--",  color = "black", label = "unit circle")
            ax.plot(X,Y, linewidth=2, label = "symbol")
            ax.axis("equal")
            ax.legend(loc="best")
            plt.title("Symbol of "+ S.name(lambda_bool = True))
        else:
            assert (type(lamb) == np.ndarray and lamb.ndim == 1) or type(lamb) == list
            plt.ylim(-1,1)
            plt.xlim(-1,1)
            for l in lamb:
                S = schem(l,order = order)
                X,Y = S.symbol(nparam)
                ax.plot(X,Y, linewidth=2, label = f"$\lambda = $ {l}")
            T = np.linspace(0,7,1000)
            ax.plot(np.cos(T),np.sin(T),"--", color = "black", label = "unit circle")
            ax.axis("equal")
            ax.legend(loc="upper left")
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
    return ax





def nbrzerosdetKL(schem, left_bound= Dirichlet(), lamb = None, sigma = 0,  order = 2, nparam=200,parametrization_bool = True, fig_size = (6,4)):
    fig,ax = plt.subplots(figsize=fig_size)
    if lamb is None:
        lambmin = 0.001
        lambmax = schem(1).CFL
    else:
        if isinstance(lamb, (int,float)):
            lamb = np.array([lamb], dtype=float)
        if len(lamb)==1:
            if lamb[0]<0:
                raise ValueError("lambda has to be non negative")
            lambmin = 0.001
            lambmax = lamb[0]
        else:
            lambmin = min(lamb)
            lambmax = max(lamb)
    lambdas = np.linspace(lambmin,lambmax,nparam)
    WindingNumbers = []
    for l in lambdas:
        S = schem(l,left_bound, sigma = sigma, order = order)
        Param, Det = S.detKL(nparam, parametrization_bool)
        WindingNumbers.append(Indice(Det))
    ax.plot(lambdas, S.r - np.array(WindingNumbers))
    plt.title(f"Number of zeros of KL determinant for {schem(1/2,left_bound, order = order, sigma=sigma).name(boundary_bool = True, sigma_bool = True, lambda_bool = False)} with respect to $\lambda$")
    return ax