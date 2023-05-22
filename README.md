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
python3 examples/draw_symbol.py
```

# Example 
```python
import matplotlib.pyplot as plt

from boundaryscheme.schemes import BeamWarming
import boundaryscheme.pyplot as bsplt

bsplt.detKLcurve(BeamWarming, SILW(2,3),lambdacursor = True)
plt.show()
```

<video src='mp4' width=180/>

# Citing

The code is citable via Zenodo. Please cite as:

P. Le Barbenchon, boundaryscheme: package Python for numerical schemes with boundaries. 2023. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7773742.svg)](https://doi.org/10.5281/zenodo.7773742)

