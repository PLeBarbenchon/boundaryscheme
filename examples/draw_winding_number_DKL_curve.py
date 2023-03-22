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
from complex_winding_number import *


#Choix du schéma 
scheme = BeamWarming
order = 5
CFL = 2
sigma = 0

#Choix du bord
boundary = SILW(2, 3)


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(6, 4))

        
    lamb = np.linspace(0.01,CFL,200)
    n_param = 60
    parametrization_bool = True

    IndiceComplexe = []
    S = scheme(1/2,boundary, sigma=sigma, order = order) 
    for l in lamb:
        S = scheme(l,boundary, sigma = sigma, order = order)
        fDKL = S.DKL()
        Param, Det = S.detKL(n_param,fDKL, parametrization_bool)
        IndiceComplexe.append(Indice(Det))

    plt.plot(lamb, S.r - np.array(IndiceComplexe))
    plt.title(f"Number of zeros of KL determinant for {scheme(1/2,boundary, order = order, sigma=sigma).name(boundary_bool = True, sigma_bool = True, lambda_bool = False)} with respect to $\lambda$")
    plt.show()

