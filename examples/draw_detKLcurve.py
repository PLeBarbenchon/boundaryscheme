"""
This file aims at testing and plotting the Kreiss-Lopatinskii determinant curve
"""

import matplotlib.pyplot as plt
import numpy as np

from boundaryscheme.boundaries import *
from boundaryscheme.schemes import *
import boundaryscheme.pyplot as bsplt


# You can import the following different schemes : BeamWarming, Upwind, LaxWendroff, LaxFriedrichs, ThirdOrder, BB, Dissipatif, LW

# You can import the following different boundaries : Dirichlet, SILW, DDJ


if __name__ == "__main__":
    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming and for lambda = 1.4"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lamb = 1.4)

    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming and for lambda in [1.4, 0.7, 1] with sigma = 0.4"""
    bsplt.detKLplot(BeamWarming, SILW(2, 3), lamb=np.array([1.4, 0.7, 1]), sigma=0)


    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming and for lambda = 1 with sigma in [-0.2, 0, 0.1,0.4]"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lamb = 1, sigma = np.array([-0.2, 0, 0.1,0.4]))

    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming and for lambda in [1.4, 0.7] with sigma in [-0.2, 0, 0.4]"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lamb = [1.4, 0.7], sigma = np.array([-0.2, 0, 0.4]))

    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming and for lambda = 1.4 and with a cursor for sigma between -0.5 and 0.5"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lamb = 1.4, sigmacursor = True)

    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming with a cursor for lambda between 0 and the CFL condition"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lambdacursor = True)

    """plot the Kreiss-Lopatinskii determinant curve for BeamWarming with a cursor for lambda between 0 and the CFL condition and a cursor for sigma between -0.5 and 0.5"""
    # bsplt.detKLplot(BeamWarming, SILW(2,3), lambdacursor = True, sigmacursor = True)
