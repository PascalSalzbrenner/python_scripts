# script to calculate the average of the quantities calculated during a ramble MD run and printed in the .track file
# written by Pascal Salzbrenner, pts28

import sys
import numpy as np

# get fileroot and the time step [ps] at which to start averaging
fileroot = sys.argv[1]
initial_step = float(sys.argv[2])

# final step for the average can also be supplied. If it is not, average until EOF
if len(sys.argv) > 3:
    # final step has been supplied
    final_step = float(sys.argv[3])
else:
    final_step = None

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

    if float(data[0]) >= initial_step and (not final_step or float(data[0]) <= final_step):
        # we are in the window where we intend to average, collect data
        volumes.append(float(data[1]))
        temperatures.append(float(data[2]))
        pressures.append(float(data[4]))
        enthalpies.append(float(data[5]))
        total_energies.append(float(data[7]))
    elif final_step and float(data[0]) > final_step:
        break
    else:
        # we haven't reached the point at which we want to average yet
        continue

trackfile.close()

if not final_step:
    final_step = float(data[0])

# calculate averages
volumes = np.array(volumes)
temperatures = np.array(temperatures)
pressures = np.array(pressures)
enthalpies = np.array(enthalpies)
total_energies = np.array(total_energies)

average_volume = volumes.mean()
average_temperature = temperatures.mean()
average_pressure = pressures.mean()
average_enthalpy = enthalpies.mean()
average_total_energy = total_energies.mean()

# calculate standard deviations
volume_std = volumes.std()
temperature_std = temperatures.std()
pressure_std = pressures.std()
enthalpy_std = enthalpies.std()
total_energy_std = total_energies.std()

# write out data
with open("{}_{:d}_{:d}_ps_averages.txt".format(fileroot, int(initial_step), int(final_step)), "w") as outfile:

    outfile.write("# Average values and their standard deviations along the MD trajectory, calculated from {} to {} ps\n\n".format(initial_step, final_step))
    outfile.write("Volume: {} A**3/atom; Standard deviation {} A**3/atom\n".format(average_volume, volume_std))
    outfile.write("Temperature: {} K; Standard deviation {} K\n".format(average_temperature, temperature_std))
    outfile.write("Pressure: {} GPa; Standard deviation {} GPa\n".format(average_pressure,  pressure_std))
    outfile.write("U_plus_PV: {} eV/atom; Standard deviation {} eV/atom\n".format(average_enthalpy, enthalpy_std))
    outfile.write("Total_Energy_Enthalpy: {} eV/atom; Standard deviation {} eV/atom\n".format(average_total_energy, total_energy_std))
