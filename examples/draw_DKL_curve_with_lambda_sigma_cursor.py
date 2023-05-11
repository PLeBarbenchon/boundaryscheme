"""
This file aims at testing and plotting the DKL function with cursor
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import os
path_here = os. getcwd()
import sys
sys.path.insert(0,path_here+"/boundaryscheme")

from schemes import *
from boundaries import *





# scheme choice
scheme = BeamWarming
order = 5
CFL = 2

# boundary choice
boundary = SILW(2, 3)


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(10, 8))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)
    ax.set_xlim(-3, 3)
    ax.axis("equal")


    
    lambda_0 = CFL/2
    sigma_0 = 0
    n_param = 100
    parametrization_bool = True

    plt.title(
        f"KL determinant for z in S for {scheme(lambda_0,boundary, order = order).shortname()} with boundary {boundary.name()}"
    )

    def calc_det(l, sigm):
        S = scheme(l, boundary, sigma=sigm, order=order)
        fDKL = S.DKL()
        return S.detKL(n_param, fDKL, parametrization_bool)[1]

    Det = calc_det(lambda_0, sigma_0)
    [line] = ax.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, color='red')
    plt.axvline(x=0, color="0.5")
    plt.axhline(y=0, color="0.5")

    axlambda = plt.axes([0.25, 0.1, 0.65, 0.03])
    lambda_slider = Slider(
        ax=axlambda,
        label='lambda',
        valmin=0.001,
        valmax=CFL,
        valinit=lambda_0,
    )

    # Make a vertically oriented slider to control the amplitude
    axsigm = plt.axes([0.1, 0.25, 0.0225, 0.63])
    sigma_slider = Slider(
        ax=axsigm, label="sigma", valmin=-1 / 2, valmax=1 / 2, valinit=sigma_0, orientation="vertical"
    )

    # The function to be called anytime a slider's value changes
    def update(val):
        Det = calc_det(lambda_slider.val, sigma_slider.val)
        line.set_xdata([z.real for z in Det])
        line.set_ydata([z.imag for z in Det])
        fig.canvas.draw_idle()

    # register the update function with each slider
    lambda_slider.on_changed(update)
    sigma_slider.on_changed(update)

    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')

    def reset(event):
        lambda_slider.reset()
        sigma_slider.reset()

    button.on_clicked(reset)

    plt.show()
