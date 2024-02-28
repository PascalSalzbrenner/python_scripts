# script to plot the phonon dispersion of a generic 3D FCC lattice
# the interactions between atoms are modelled by harmonic springs connecting nearest neighbours
# this is analogous to a first nearest-neighbour Born-von Karman model
# the eigenvalue problem and its solution are taken from https://lampz.tugraz.at/~hadley/ss1/phonons/fcc/fcc.php

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
import matplotlib.pyplot as plt

from numpy.polynomial import Polynomial

# define the function returning the third order polynomial

def poly_equation(kx, ky, kz):
    m11 = 4 - np.cos(kx/2 + ky/2) - np.cos(kx/2 + kz/2) - np.cos(kx/2 - ky/2) - np.cos(kx/2 - kz/2)
    m12 = - np.cos(kx/2 + ky/2) + np.cos(kx/2 - ky/2)
    m13 = - np.cos(kx/2 + kz/2) + np.cos(kx/2 - kz/2)
    m22 = 4 - np.cos(kx/2 +ky/2) - np.cos(ky/2 + kz/2) - np.cos(kx/2 - ky/2) - np.cos(ky/2 - kz/2)
    m23 = - np.cos(ky/2 + kz/2) + np.cos(ky/2 - kz/2)
    m33 = 4 - np.cos(kx/2 + kz/2) - np.cos(ky/2 + kz/2) - np.cos(kx/2 - kz/2) - np.cos(ky/2 - kz/2)
    
    b = m11 + m22 + m33
    c = m12*m12 + m23*m23 - m11*m22 - m11*m33 - m22*m33
    d = m11*m22*m33 + 2*m12*m13*m23 - m12*m12*m33 - m13*m13*m22 - m23*m23*m11

    return Polynomial([d, c, b, -1])

# define dictionary of high-symmetry k-points
hsp = {"G": np.array([0, 0, 0]), "X": np.array([2*np.pi, 0, 0]), "L": np.array([np.pi, np.pi, np.pi]), "W": np.array([2*np.pi, np.pi, 0]), "U": np.array([2*np.pi, np.pi, np.pi]), "K": np.array([3*np.pi/2, 3*np.pi/2, 0])}

# define the path in k-space
k_path = ["G", "X", "W", "K", "G", "L"]

# set up output lists - linear k-path on the x-axis, 3 frequencies on the y-axis
k_linear = []
omega1 = []
omega2 = []
omega3 = []

# define list for x-axis labels
hsp_linear = []

# set initial offset
offset = 0

# loop over the k-path
for i in range(len(k_path)-1):

    # calculat the k-point coordinates
    k1 = hsp[k_path[i]]
    k2 = hsp[k_path[i+1]]

    # calculat the distance between the k-points
    dk = np.linalg.norm(k2 - k1)

    # define the number of points between the k-points
    n = 10000

    # loop over the points between the k-points
    for j in range(n):

        # calculate the k-point
        k = k1 + (k2 - k1) * j / n

        # construct the third order polynomial
        c = poly_equation(k[0], k[1], k[2])

        # get the roots of the polynomial
        r = c.roots()

        # sort the roots
        r = np.sort(r)

        # calculate offset to the k-point and add labels
        if j == n-1:
            hsp_linear.append([offset, k_path[i]])
            offset += dk
            if i == len(k_path)-2:
                hsp_linear.append([offset, k_path[i+1]])

        # append the k-point to the list
        k_linear.append(np.linalg.norm(k - k1) + offset)

        # append the frequencies to the list
        omega1.append(r[0])
        omega2.append(r[1])
        omega3.append(r[2])

# plot the phonon dispersion
colour = "#8000C4"
plt.plot(k_linear, omega1, color=colour)
plt.plot(k_linear, omega2, color=colour)
plt.plot(k_linear, omega3, color=colour)

# labels etc

# Hide x-axis and y-axis
plt.xticks([])

# y-axis label
plt.ylabel(r'$\sqrt{\frac{M}{C}}\omega$', fontsize=14)

# x-axis label and bars
for point in hsp_linear:
    x = point[0]
    y = point[1]
    plt.axvline(x=x, color='black', linestyle=':')
    plt.text(x, -0.1, f'{y}', color='black', ha='center')

# save
plt.savefig("fcc_phonon_dispersion.png", dpi=300)
