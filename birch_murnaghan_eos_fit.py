# short script to fit to the Birch-Murnaghan equation of state (DOI: 10.1103/PhysRev.71.809)
# input: energy.txt file as produced by my geometry_optimisation_cubic routine

# written by Pascal Salzbrenner - pts28@cam.ac.uk

import sys
import math
import numpy as np
from scipy.optimize import curve_fit

# define Birch-Murnaghan equation of state for E(V) - V is input, the other variables are
# parameters to be fit
def birch_murnaghan_eos(V, E_0, V_0, B_0, B_prime):
    """The Birch-Murnaghan equation of state as defined in DOI: 10.1103/PhysRev.71.809"""

    return E_0 + 9.0*V_0*B_0/(16.0) * ((((V_0/V)**(2))**(1.0/3.0)-1)**(3)*B_prime
           + (((V_0/V)**(2))**(1.0/3.0) - 1)**(2)*(6.0-4.0*((V_0/V)**(2))**(1.0/3.0)))

# determine the prefactor
spacegroup = sys.argv[1]

if spacegroup.startswith("P"):
    # simple cubic
    prefactor = 1
elif spacegroup.startswith("F"):
    # FCC
    prefactor = 0.5
elif spacegroup.startswith("I"):
    # BCC
    prefactor = 0.25
else:
    prefactor = 1

# conversion factor from Angstrom to Bohr
angstrom_to_bohr = 1.8897261254578281

# define x and y data
volumes = []
energies = []

# read energy and volume out of the energy.txt file
with open("energy.txt", "r") as energyfile:

    for line in energyfile:
        if line.lstrip().startswith("#"):
            # comment line
            continue
        volumes.append(float(line.split()[1])) # volume
        energies.append(float(line.split()[2])) # energy

parameters, covariance = curve_fit(birch_murnaghan_eos, volumes, energies,
bounds=((2*min(energies), volumes[0], 0, -np.inf), (max(energies), volumes[-1], np.inf, np.inf)),
maxfev=100000) # impose very general boundaries to ensure correct result
# if the sampling in geometry_optimisation_cubic is anywhere near reasonable,
# this reliably gives the right result

with open("birch_murnaghan_eos_parameters.txt", "w") as outfile:

    outfile.write("E_0 = {} eV\n".format(parameters[0]))
    outfile.write("V_0 = {} A^3 = {} Bohr^3\n".format(parameters[1], parameters[1]*angstrom_to_bohr**3))
    outfile.write("lattice parameter a = {} A = {} Bohr\n".format((prefactor*np.abs(parameters[1]))**(1.0/3.0), (prefactor*np.abs(parameters[1]))**(1.0/3.0)*angstrom_to_bohr))
    outfile.write("B_0 = {} eV/A^3 = {} GPa\n".format(parameters[2], parameters[2]*1.60217662*10*10)) # conversion eV/A^3 to GPa
    outfile.write("B_prime = {} [dimensionless]".format(parameters[3]))

# add writing of gnuplot file to plot fit
with open("energy_volume.gnu", "w") as plotfile:

    plotfile.write("set terminal postscript eps colour font 'Helvetica,20'\n")
    plotfile.write("set style data points\n")
    plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque font 'Helvetica,15'\n")
    plotfile.write("set output '| epstopdf --filter --outfile=energy_volume.pdf'\n")
    plotfile.write("E_0 = {}\n".format(parameters[0]))
    plotfile.write("V_0 = {}\n".format(parameters[1]))
    plotfile.write("B_0 = {}\n".format(parameters[2]))
    plotfile.write("B_prime = {}\n".format(parameters[3]))
    plotfile.write("birch_murnaghan_eos(x) = E_0 + 9.0*V_0*B_0/(16.0)*(((V_0/x)**(2.0/3.0)-1)**(3)*B_prime + ((V_0/x)**(2.0/3.0)-1)**(2)*(6-4*(V_0/x)**(2.0/3.0)))\n")
    plotfile.write("set xrange [{}:{}]\n".format(math.floor(min(volumes)), math.ceil(max(volumes))))
    plotfile.write("set xlabel 'Volume [A^3]'\n")
    plotfile.write("set ylabel 'Energy [eV]'\n")
    plotfile.write("plot 'energy.txt' u 2:3 w points pt 7 ps 1.5 lc rgb '#DC143C' title 'Energy-Volume data from DFT calculation', birch_murnaghan_eos(x) lw 2 lc rgb '#DC143C' title 'Least-Squares fit to Birch-Murnaghan EOS'")
