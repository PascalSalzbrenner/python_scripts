# script to plot the phonon dispersion of a generic 3D FCC lattice
# the interactions between atoms are modelled by harmonic springs connecting nearest neighbours
# this is analogous to a first nearest-neighbour Born-von Karman model
# the eigenvalue problem and its solution are taken from https://lampz.tugraz.at/~hadley/ss1/phonons/fcc/fcc.php

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
import matplotlib.pyplot as plt

from numpy.polynomial import Polynomial

# define the function returning the coefficient of the third order polynomial

