"""
This file aims at testing and plotting the DKL function
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import boundaryscheme as bs
from boundaryscheme.schemes import BeamWarming
import boundaryscheme.pyplot as pl


if __name__ == "__main__":
    pl.symbolplot(BeamWarming, lambdacursor = True)

    # bs.symbolplot(BeamWarming, [0.7,1.4, 1.8], lambdacursor = True)

    # bs.symbolplot(BeamWarming, 1.4)

    # bs.symbolplot(BeamWarming, np.array([0.7, 1, 1.6]))
