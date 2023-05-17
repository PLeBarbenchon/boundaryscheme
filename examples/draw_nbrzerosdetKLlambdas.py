"""
This file aims at testing and plotting the number of zeros of the Kreiss-Lopatinskii determinant with respect to  lambda
""" 

import matplotlib.pyplot as plt
import numpy as np

from boundaryscheme.schemes import BeamWarming
from boundaryscheme.boundaries import SILW
import boundaryscheme.pyplot as bsplt


#You can import the following different schemes : BeamWarming, Upwind, LaxWendroff, LaxFriedrichs, ThirdOrder, BB, Dissipatif, LW

#You can import the following different boundaries : Dirichlet, SILW, DDJ


if __name__ == "__main__":
    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with SILW(2,3) for lambda between 0 and the CFL condition"""
    bsplt.nbrzerosdetKL(BeamWarming, SILW(2,3))

    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with SILW(2,3) for lambda between 1.4 and 1.7"""
    # bsplt.nbrzerosdetKL(BeamWarming, SILW(2,3),lamb = [1.4,1.7], nparam = 100)

    plt.show()

