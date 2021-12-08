#postprocessing utility to easily read the highest/lowest value of a specified band from Quantum Espresso bands.x output
#written by Pascal Salzbrenner - pts28@cam.ac.uk

import os
import numpy as np

filename = input("Please enter the name of the file containing the band structure data: ")
band_index = int(input("Please enter the index of the band whose maximum and minimum you want to evaluate: "))

band_data = open("{}".format(filename), "r")

#read past the header line
line = band_data.readline()

# all following lines give alternatingly the k-point and then the energy eigenvalues at it

# treat the first set uniquely to set the baselines

k_vector = np.array(band_data.readline().split(), dtype=float)
energies = np.array(band_data.readline().split(), dtype=float)

#the k-points associated with the highest and lowest energies of the band
k_point_E_max = k_vector[0:3]
k_point_E_min = k_vector[0:3]

E_max = energies[band_index-1]
E_min = energies[band_index-1]

for k_point in band_data:
    k_vector = np.array(k_point.split(), dtype=float)
    energies = np.array(band_data.readline().split(), dtype=float)

    if float(energies[band_index-1]) > E_max:
        E_max = float(energies[band_index-1])
        k_point_E_max = k_vector[0:3]
    elif float(energies[band_index-1]) < E_min:
        E_min = float(energies[band_index-1])
        k_point_E_min = k_vector[0:3]

band_data.close()

#write output
with open("band_{}_max_min.dat".format(str(band_index)), "w") as outfile:
    outfile.write("The band's highest energy is {} eV. In units of the reciprocal lattice vectors, it is at: {} {} {}\n".format(E_max, k_point_E_max[0],
    k_point_E_max[1], k_point_E_max[2]))
    outfile.write("The band's lowest energy is {} eV. In units of the reciprocal lattice vectors, it is at: {} {} {}\n".format(E_min, k_point_E_min[0],
    k_point_E_min[1], k_point_E_min[2]))
