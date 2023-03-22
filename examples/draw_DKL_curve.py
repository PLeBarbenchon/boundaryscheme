"""
This file aims at testing and plotting the DKL function
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# appending a path
sys.path.append('thesis')
from schemes import *
from boundaries import *


#Choix du schéma 
scheme = BeamWarming
order = 5
CFL = 2
sigma = 0
lamb = np.array([0.3,1.4,1.95])

#Choix du bord
boundary = SILW(2, 3)


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(6, 4))


    n_param = 300
    parametrization_bool = False

    for l in lamb:
        S = scheme(l, boundary, order=order, sigma=sigma)
        fDKL = S.DKL()
        Param, Det = S.detKL(n_param, fDKL, parametrization_bool)
        plt.plot([z.real for z in Det], [z.imag for z in Det], linewidth=2, label=f"$\lambda$ = " + str(l))

    plt.axvline(x=0, color="0.5")
    plt.axhline(y=0, color="0.5")
    plt.legend(loc='best')

    if sigma != 0:
        plt.title(
            f"KL determinant for z in S for {S.shortname()} with boundary {boundary.name()} with $\sigma = ${sigma}"
        )
    else:
        plt.title(f"KL determinant for z in S for {S.shortname()} with boundary {boundary.name()}")

    # plt.xlim(-9*10**(-6),9*10**(-6))
    # plt.ylim(-6*10**(-6),6*10**(-6))
    plt.axis('equal')
    plt.show()
