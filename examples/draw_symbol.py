"""
This file aims at testing and plotting the DKL function
"""

import numpy as np
import matplotlib.pyplot as plt
import os
path_here = os. getcwd()
import sys
sys.path.insert(0,path_here+"/boundaryscheme")

from schemes import *
import pyplot as bs



if __name__ == "__main__":
    bs.symbolplot(BeamWarming, lambdacursor = True)

    bs.symbolplot(BeamWarming, [0.7,1.4, 1.8], lambdacursor = True)

    bs.symbolplot(BeamWarming, 1.4)

    bs.symbolplot(BeamWarming, np.array([0.7, 1, 1.6]))
