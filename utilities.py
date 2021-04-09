# Utilities which are useful for several of my scripts

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

########################################### conversion between Cartesian and direct coordinates ###########################################

def lattice_basis_to_cartesian(point, lattice):
    """Function to convert an atom's coordinates, given in direct coordinates, into Cartesian coordinates.
    Concretely, with the units of the lattice
    :param np.ndarray point: numpy array of length dimensions indicating a point in the lattice basis
    :param np.ndarray lattice: numpy array of size dimensions x dimensions, containing the lattice vectors as rows

    :returns np.ndarray point_cart: numpy array of length dimensions giving the point in Cartesian coordinates
    """

    return np.dot(lattice.T, point)

def cartesian_to_lattice_basis(point, lattice):
    """Function to convert a point given in Cartesian coordinates into coordinates in the lattice
    :param np.ndarray point: numpy array of length dimensions indicating a point in Cartesian coordinates
    :param np.ndarray lattice: numpy array of size dimensions x dimensions, containing the lattice vectors as rows

    :returns np.ndarray point_lattice_basis: numpy array of length dimensions giving the point in the lattice basis
    """

    return np.dot(np.linalg.inv(lattice.T), point)
