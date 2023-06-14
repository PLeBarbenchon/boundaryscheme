This is the architecture for the package "boundaryscheme" : https://pypi.org/project/boundaryscheme/


# boundaryscheme
Package Python to use numerical scheme with boundaries which is described in the PhD manuscript 
> P. Le Barbenchon, Étude théorique et numérique de la stabilité GKS pour des schémas d'ordre élevé en présence de bords, PhD, 2023.

# Package documentation

https://plebarbenchon.github.io/boundaryscheme

# Easy installation for the PyPI version

```bash
pip install boundaryscheme
```

# Installation for the GitHub version
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

# Creation of a scheme

Pour écrire le schéma $U_{j}^{n+1} = \dfrac{1}{3}U_{j-1}^n + U_j^n -\dfrac{1}{6} U_{j+1}^n - \dfrac{1}{6}U_{j+2}^n$ avec un certain bord $B$, il faut écrire 

```python
S = Scheme([1/3,1,-1/6,-1/6],1, B) 
```


Pour écrire un schéma qui dépend d'un paramètre $\lambda$, on peut créer la classe du schéma de la façon suivante :

```python
class Name(Scheme):
    """This is a class to represent ...

    :param lamb: The Courant number, i.e  a.dt/dx where "a" is the velocity, "dt" the time discretization and "dx" the space discretization
    :type lamb: float
    :param boundary: Boundary condition, defaults to Dirichlet()
    :type boundary: class:`Boundary`, optional
    :param sigma: Gap between the mesh and the boundary condition, defaults to 0
    :type sigma: float, optional
    """

    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, D=0, **kwargs):
        """Constructor method"""
        self.sigma = sigma
        self.lamb = lamb
        self.inter = #write the list of the coefficients of the scheme
        self.center = #write the index of the center of the scheme
        self.CFL = #give the CFL condition
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "Name"
```
Par exemple pour écrire le schéma 
$U_{j}^{n+1} = \dfrac{\lambda}{3}U_{j-1}^n + \lambda^2 U_j^n -\dfrac{1}{6} U_{j+1}^n - \dfrac{\lambda}{6}U_{j+2}^n$ avec un certain bord $B$, il faut écrire 

```python
class Name(Scheme):
    def __init__(self, lamb, boundary=Dirichlet(), sigma=0, D=0, **kwargs):
        self.sigma = sigma
        self.lamb = lamb
        self.inter = [lamb/3, lamb**2, -1/6, -lamb/6]
        self.center = 1
        self.CFL = #give the CFL condition
        super().__init__(inter=self.inter, center=self.center, boundary=boundary, sigma=sigma, **kwargs)

    def shortname(self):
        """Name method"""
        return "Name"
```


# Citing

The code is citable via Zenodo. Please cite as:

P. Le Barbenchon, boundaryscheme: package Python for numerical schemes with boundaries. 2023. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7773742.svg)](https://doi.org/10.5281/zenodo.7773742)

