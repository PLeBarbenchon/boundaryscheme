"""
This file aims at testing and plotting the number of zeros of the Kreiss-Lopatinskii determinant
""" 

import matplotlib.pyplot as plt
import numpy as np

from boundaryscheme.schemes import *
from boundaryscheme.boundaries import *
import boundaryscheme.pyplot as bsplt


#You can import the following different schemes : BeamWarming, Upwind, LaxWendroff, LaxFriedrichs, ThirdOrder, BB, Dissipatif, LW

#You can import the following different boundaries : Dirichlet, SILW, DDJ


if __name__ == "__main__":
    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with SILW(2,3) for lambda between 0 and the CFL condition"""
    # bsplt.nbrzerosdetKL(BeamWarming, SILW(2,3))

    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with SILW(2,3) for lambda between 1.4 and 1.7"""
    # bsplt.nbrzerosdetKL(BeamWarming,  SILW(2,3),lamb = np.array([1.4,1.7]), nparam = 150)


    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with several boundary conditions"""
    # bsplt.nbrzerosdetKL(BeamWarming, [SILW(1,3), SILW(2,3),SILW(2,4),SILW(3,4),DDJ(3,0),DDJ(3,1)],nlambda = 60)

    """plot the number of zeros of Kreiss-Lopatinskii determinant of BeamWarming with SILW(2,3) with respect to lambda and sigma"""
    bsplt.nbrzerosdetKL(BeamWarming, SILW(2,3), sigma=True, nlambda=50, nsigma = 50)

    plt.show()
 

