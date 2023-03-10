"""
This file aims at testing and plotting the DKL function
"""

import numpy as np
import matplotlib.pyplot as plt

from thesis.schemes import LW
from thesis.boundaries import DDJ


if __name__ == "__main__":
    fig = plt.figure(1, figsize=(6, 4))

    scheme = LW
    boundary = DDJ(6, 1)

    order = 5
    sigma = 0
    lamb = np.array([0.01])
    n_param = 300
    parametrization_bool = False

    for l in lamb:
        S = scheme(l, boundary, order=order, sigma=sigma)
        fDKL = S.DKL()
        Det = S.detKL(n_param, fDKL, parametrization_bool)
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
