This is the architecture for the package "boundaryscheme" 


# boundaryscheme
Package Python to use numerical scheme with boundaries

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
from boundaryscheme.schemes import BeamWarming
import boundaryscheme.pyplot as bsplt

bsplt.symbolplot(BeamWarming, lambdacursor = True)
plt.show()
```

# Citing

The code is citable via Zenodo. Please cite as:

P. Le Barbenchon, boundaryscheme: package Python for numerical schemes with boundaries. 2023. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7773742.svg)](https://doi.org/10.5281/zenodo.7773742)
