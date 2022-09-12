# script to calculate the average of the quantities calculated during a ramble MD run and printed in the .track file
# written by Pascal Salzbrenner, pts28

import sys
import numpy as np

# get fileroot and the time step [ps] at which to start averaging
fileroot = sys.argv[1]
initial_step = float(sys.argv[2])

# open track file
trackfile = open("{}.track".format(fileroot), "r")

# set up lists for the different quantities
volumes = [] # [A**3/atom], second column
temperatures = [] # [K], third column
pressures = [] # [GPa], fifth column
enthalpies = [] # [eV/atom], sixth column
total_energies = [] # (enthalpy + kinetic energy) [eV/atom], eighth column

# read data
for line in trackfile:
    data = line.split()

    if float(data[0]) >= initial_step:
        # we are at or beyond the time step from which we intend to average, collect data
        volumes.append(float(data[1]))
        temperatures.append(float(data[2]))
        pressures.append(float(data[4]))
        enthalpies.append(float(data[5]))
        total_energies.append(float(data[7]))
    else:
        # we haven't reached the point at which we want to average yet
        continue

trackfile.close()

# calculate averages
volumes = np.array(volumes)
temperatures = np.array(temperatures)
pressures = np.array(pressures)
enthalpies = np.array(enthalpies)
total_energies = np.array(total_energies)

average_volume = volumes.mean
average_temperature = temperatures.mean
average_pressure = pressures.mean
average_enthalpy = enthalpies.mean
average_total_energy = total_energies.mean

# write out data
with open("{}_averages.txt".format(fileroot), "w") as outfile:

    outfile.write("# Average values along the MD trajectory, calculated from {} to {} ps\n\n".format(initial_step, data[0]))
    outfile.write("Volume: {} A**3/atom\n".format(average_volume))
    outfile.write("Temperature: {} K\n".format(average_temperature))
    outfile.write("Pressure: {} GPa\n".format(average_pressure))
    outfile.write("Enthalpy: {} eV/atom\n".format(average_enthalpy))
    outfile.write("Total_Energy: {} eV/atom\n".format(average_total_energy))
