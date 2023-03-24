"""
This file aims at testing and plotting the DKL function with cursor
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import os
path_here = os. getcwd()
import sys
sys.path.insert(0,path_here+"/thesis")

from schemes import *
from boundaries import *

# scheme choice
scheme = BeamWarming
order = 5
CFL = 2
sigma = 0

# boundary choice
boundary = SILW(2, 3)

if __name__ == "__main__":
    fig = plt.figure(1, figsize=(10, 8))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    lambda_0 = CFL / 2
    n_param = 100
    parametrization_bool = True

    if sigma != 0:
        plt.title(
            f"KL determinant for z in S for {scheme(lambda_0,boundary, order = order).shortname()} with boundary {boundary.name()} with $\sigma = ${sigma}"
        )
    else:
        plt.title(
            f"KL determinant for z in S for {scheme(lambda_0,boundary, order = order).shortname()} with boundary {boundary.name()}"
        )

    def calc_det(l):
        S = scheme(l, boundary, sigma=sigma, order=order)
        fDKL = S.DKL()
        return S.detKL(n_param, fDKL, parametrization_bool)[1]

    Det = calc_det(lambda_0)
    [line] = ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, color='red')
    plt.axvline(x=0, color="0.5")
    plt.axhline(y=0, color="0.5")
    ax.axis("equal")

    lambda_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
    lambda_slider = Slider(lambda_slider_ax, 'lambda', 0.0001, CFL, valinit=lambda_0)

    def sliders_on_changed(val):
        Det = calc_det(lambda_slider.val)
        line.set_xdata([z.real for z in Det])
        line.set_ydata([z.imag for z in Det])
        fig.canvas.draw_idle()

    lambda_slider.on_changed(sliders_on_changed)

    plt.show()
