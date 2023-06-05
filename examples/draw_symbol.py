"""
This file aims at testing and plotting the symbol curve
"""

import matplotlib.pyplot as plt
import numpy as np

from boundaryscheme.schemes import *
import boundaryscheme.pyplot as bsplt


# You can import the following different schemes : BeamWarming, Upwind, LaxWendroff, LaxFriedrichs, ThirdOrder, BB, Dissipatif, LW


if __name__ == "__main__":
    """plot the BeamWarming symbol curve with a cursor for lambda between 0 and the CFL condition"""
    bsplt.symbolplot(BeamWarming, lambdacursor=True)

    """plot the BeamWarming symbol curve with a cursor for lambda between 0.7 and 1.8"""
    # bsplt.symbolplot(BeamWarming, [0.7,1.4, 1.8], lambdacursor = True)

    """plot the BeamWarming symbol curve for lambda = 1.4"""
    # bsplt.symbolplot(BeamWarming, 1.4)

    """plot the BeamWarming symbol curves for lambda in [0.7, 1, 1.6]"""
    # bsplt.symbolplot(BeamWarming, np.array([0.7, 1, 1.6]))

    plt.show()
