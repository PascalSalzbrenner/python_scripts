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
    c = m12*m12 + m13*m13 + m23*m23 - m11*m22 - m11*m33 - m22*m33
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

# loop over the k-path
for i in range(len(k_path)-1):

    # set offset
    if i == 0:
        offset = 0
    else:
        offset += dk

    # calculat the k-point coordinates
    k1 = hsp[k_path[i]]
    k2 = hsp[k_path[i+1]]

    # calculat the distance between the k-points
    dk = np.linalg.norm(k2 - k1)

    # add to labels list
    hsp_linear.append([offset, k_path[i]])

    # define the number of points between the k-points
    n = 1000

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

        # append the k-point to the list
        k_linear.append(np.linalg.norm(k - k1) + offset)

        # append the frequencies to the list
        omega1.append(np.sqrt(r[0]))
        omega2.append(np.sqrt(r[1]))
        omega3.append(np.sqrt(r[2]))

    if i == len(k_path)-2:
        hsp_linear.append([offset+dk, k_path[i+1]])

        k=k2


        # construct the third order polynomial
        c = poly_equation(k[0], k[1], k[2])

        # get the roots of the polynomial
        r = c.roots()

        # sort the roots
        r = np.sort(r)       

        # append the k-point to the list
        k_linear.append(np.linalg.norm(k - k1) + offset)

        # append the frequencies to the list
        omega1.append(np.sqrt(r[0]))
        omega2.append(np.sqrt(r[1]))
        omega3.append(np.sqrt(r[2]))

# plot the phonon dispersion
colour = "#8000C4"
plt.plot(k_linear, omega1, color=colour)
plt.plot(k_linear, omega2, color=colour)
plt.plot(k_linear, omega3, color=colour)

# labels etc

# Hide x-axis
plt.xticks([])

plt.xlim(hsp_linear[0][0], hsp_linear[-1][0])
plt.ylim(0, 3)

# y-axis label
plt.ylabel(r'$\sqrt{\frac{M}{C}}\omega$', fontsize=15)

# x-axis label and bars
for point in range(len(hsp_linear)):
    x = hsp_linear[point][0]
    name = hsp_linear[point][1]
    if name == "G":
        plt.text(x, -0.2, "Γ", color='black', ha='center', fontsize=15)
    else:
        plt.text(x, -0.2, f'{name}', color='black', ha='center', fontsize=15)

    if point != 0 and point != len(hsp_linear)-1:
        plt.axvline(x=x, color='black', linestyle='-')

# save
plt.savefig("fcc_phonon_dispersion.png", dpi=300)
