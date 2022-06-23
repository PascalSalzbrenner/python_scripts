# script to fit a parabola to the energy as a function of the ratio x (number of guest atoms per 8 host atoms) for host-guest structures
# for the definition, see the supplementary material of DOI: 10.1038/nmat2796
# inputs: seed.res files for each structure, which contain the energies;
# seed_x.dat files corresponding to the .res files - these contain the x values

# written by Pascal Salzbrenner - pts28@cam.ac.uk

import os
import sys
import numpy as np
from numpy.polynomial import Polynomial

units = sys.argv[1].lower() # eV or meV

if len(sys.argv) == 1 or sys.argv[1].lower() != "mev":
    # unless meV is explicitly specified, we use eV
    factor = 1
    units = "eV"
else:
    factor = 1000
    units = "meV"

# detect all structures for which a seed_x.dat file is present and read the data
x_list = []
energies = [] # in eV/atom

for file in os.listdir():
    if "_x.dat" in file:
        # read value of x
        with open(file, "r") as x_file:
            x_list.append(float(x_file.readline()))

        # read energy and n_atoms and calculate energy per atom
        with open("{}.res".format(file.rstrip("_x.dat")), "r") as resfile:
            title_line = resfile.readline().split()
            energies.append(float(title_line[4])/float(title_line[7]))

x_list = np.array(x_list)
energies = np.array(energies)

# fit the parabola
quadratic_fit = Polynomial.fit(x_list, energies, 2, window=[x_list.min(), x_list.max()])
coefficients = quadratic_fit.coef

# find minimum - where the first derivative is 0
quadratic_fit_derivative = quadratic_fit.deriv(1)
min = quadratic_fit_derivative.roots()

# write the minimum to a file
with open("optimal_x.txt", "w") as x_file:
    x_file.write("#x which minimises the energy; energy [eV/atom]\n")
    x_file.write("{} {}\n".format(min[0], min[1]))

# write x - energy data to a file
with open("x_energy.dat", "w") as energy_file:
    energy_file.write("# x; energy [{}/atom]\n".format(units))
    for i in range(len(x_list)):
        energy_file.write("{} {}\n".format(x_list[i], energies[i]*factor-min[1]*factor))

# plot data
plotfile = open("x_energy.gnu", "w") # gnuplot file to plot fit

plotfile.write("set terminal postscript eps colour font 'Helvetica,20'\n")
plotfile.write("set style data points\n")
plotfile.write("set output '| epstopdf --filter --outfile=x_energy.pdf'\n")

plotfile.write("set xlabel 'x'\n")
plotfile.write("set ylabel 'Enthalpy [{}/atom]'\n".format(units))
plotfile.write("set xrange [{}:{}]\n".format(x_list.min(), x_list.max()))

plotfile.write("e(x)={}*x**2 + {}*x + {}".format(coefficients[2], coefficients[1], coefficients[0]))

plotfile.write("plot e(x) lc rgb '#000080', 'x_energy.dat' u 1:2 pt 7 lc rgb '#DC143C', 'optimal_x.txt' u 1:2 pt 7 lc rgb 'black'\n")

plotfile.close()
