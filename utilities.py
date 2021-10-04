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

################ generation of a unit cell from the abcABC format (lengths of lattice vectors and the angles between them) ################

def construct_lattice_from_abc(vector_lengths, angles):
    """Given lattice data in the a b c alpha beta gamma format, this function constructs the lattice vectors
       see the explanation for determining lattice vectors from this format at
       https://en.wikipedia.org/wiki/Fractional_coordinates#In_crystallography

       :param list vector_lengths: a list containing [a, b, c]
       :param list angles: a list containing [alpha, beta, gamma] in degrees
    """

    v_1 = [float(vector_lengths[0]), 0, 0]
    v_2 = [float(vector_lengths[1])*np.cos(float(angles[2])*degree_to_radian),
    float(vector_lengths[1])*np.sin(float(angles[2])*degree_to_radian), 0]
    v_3_x = float(vector_lengths[2])*np.cos(float(angles[1])*degree_to_radian)
    v_3_y = float(vector_lengths[2])*(np.cos(float(angles[0])*degree_to_radian)-np.cos(float(angles[2])*degree_to_radian)
            *np.cos(float(angles[1])*degree_to_radian))/(np.sin(float(angles[2])*degree_to_radian))
    v_3_z = np.sqrt(float(vector_lengths[2])**2-v_3_x**2-v_3_y**2)

    lattice = np.array([v_1, v_2, [v_3_x, v_3_y, v_3_z]])

    return lattice
