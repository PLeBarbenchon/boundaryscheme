"""
This file aims at testing and plotting the DKL function
"""

from boundaryscheme.boundaries import SILW
from boundaryscheme.schemes import BeamWarming
import boundaryscheme.pyplot as bsplt


if __name__ == "__main__":
    bsplt.symbolplot(BeamWarming, lambdacursor = True)
    bsplt.show()

    # bsplt.symbolplot(BeamWarming, [0.7,1.4, 1.8], lambdacursor = True)

    # bsplt.symbolplot(BeamWarming, 1.4)

    # bsplt.symbolplot(BeamWarming, np.array([0.7, 1, 1.6]))
