This is the architecture for the package "boundaryscheme" 


# boundaryscheme
Package Python to use numerical scheme with boundaries which is described in the PhD manuscript 
> P. Le Barbenchon, Étude théorique et numérique de la stabilité GKS pour des schémas d'ordre élevé en présence de bords, PhD, 2023.

# Installation
```bash
git clone https://github.com/PLeBarbenchon/boundaryscheme.git
cd boundaryscheme
pip3 install -r requirements.txt
pip3 install -e .
python3 examples/draw_detKLcurve.py
```

# Example 
```python
import matplotlib.pyplot as plt

from boundaryscheme.schemes import BeamWarming
from boundaryscheme.boundaries import SILW
import boundaryscheme.pyplot as bsplt

bsplt.detKLcurve(BeamWarming, SILW(2,3),lambdacursor = True)
plt.show()
```

![mygif](https://github.com/PLeBarbenchon/boundaryscheme/assets/92107096/2ca0d414-77a6-410e-a582-a3950699dcf0)

# Plot a Kreiss--Lopatinskii determinant curve

à écrire

# Create a scheme

In schemes.py, create 

# Citing

The code is citable via Zenodo. Please cite as:

P. Le Barbenchon, boundaryscheme: package Python for numerical schemes with boundaries. 2023. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7773742.svg)](https://doi.org/10.5281/zenodo.7773742)

