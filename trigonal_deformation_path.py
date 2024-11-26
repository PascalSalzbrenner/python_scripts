# script to interpolate the singlepoint energies calculated along the trigonal deformation path
# (fcc -> R-3m -> sc -> R-3m -> bcc) for a series of angles and interatomic distances

# finds the optimum interatomic distance for a series of pressures and angles
# this is done using the Birch-Murnaghan equation of state, DOI: 10.1103/PhysRev.71.809

# the output .res files have energies only, so we add the PV term manually for the pressure range

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


################################################ Definitions and setup ################################################

# define Birch-Murnaghan equation of state for E(V) - V is input, the other variables are
# parameters to be fit
def birch_murnaghan_eos(V, E_0, V_0, B_0, B_prime):
    """The Birch-Murnaghan equation of state as defined in DOI: 10.1103/PhysRev.71.809"""

    return E_0 + 9.0*V_0*B_0/(16.0) * ((((V_0/V)**(2))**(1.0/3.0)-1)**(3)*B_prime
           + (((V_0/V)**(2))**(1.0/3.0) - 1)**(2)*(6.0-4.0*((V_0/V)**(2))**(1.0/3.0)))

# define pressure conversion factor from eV/A**3 to GPa
# 1 eV/A**3 = 160.21766208 GPa
pressure_conversion = 160.21766208

# get pressure range in GPa
min_press = int(sys.argv[1])
max_press = int(sys.argv[2])
press_step = int(sys.argv[3])

# in the end, we need lists of angles, optimum volumes at the respective angles, and corresponding optimum energies
# once for each pressure, so we set up a dictionary for this
# entries are pressure: [angles], [opt_volumes], [opt_energies]
output_dict = {}

# not elegant, but we also need a pressure list for plotting
press_list = []

# dictionary for all the energies and volumes at a given angle
data_dict = {}

################################################### Read input data ###################################################

ls = sorted([file for file in os.listdir() if file.endswith(".res")])

# populate dictionary with data
for file in ls:

    # note that the angles here will be str, which is useful for dict keys
    angle = file.split("-")[1]

    if not round(float(angle), 5) in data_dict.keys():
        # first list volumes, second list single-point energies
        data_dict[round(float(angle), 5)] = [[], []]

    with open(file) as f:
        data = f.readline().split()
        data_dict[round(float(angle), 5)][0].append(float(data[3]))
        data_dict[round(float(angle), 5)][1].append(float(data[4]))

################################## Minimisation of energy for each pressure and angle ##################################

for press in range(min_press, max_press+1, press_step):

    press_list.append(press)

    output_dict[press] = [[], [], []]

    for angle, data in sorted(data_dict.items()):
        volumes = np.array(data[0])
        energies = (np.array(data[1]) + volumes * press/pressure_conversion)#*1000
        #energies -= np.min(energies)
        parameters, covariance = curve_fit(birch_murnaghan_eos, volumes, energies,
                                           bounds=((2 * min(energies), volumes[0], 0, -np.inf),
                                                   (max(energies), volumes[-1], np.inf, np.inf)), maxfev=100000)

        output_dict[press][0].append(angle)
        output_dict[press][1].append(parameters[1])
        output_dict[press][2].append(parameters[0])

        # plot spot checks

        if int(angle) % 10 == 0:
            fit_volumes = np.arange(volumes[0], volumes[-1], 0.001)
            fit_energies = []

            for fit_volume in fit_volumes:
                fit_energies.append(birch_murnaghan_eos(fit_volume, parameters[0], parameters[1], parameters[2], parameters[3]))

            plt.plot(volumes, energies, 'o', fit_volumes, fit_energies, '-')
            plt.savefig("energy_volume_birch_murnaghan_fit_{}_GPa_{}_degree.pdf".format(press, angle))
            plt.close()

# optionally, could interpolate the minimum energy as a function of angle

####################################################### plotting #######################################################

fig, axes = plt.subplots(nrows=len(output_dict.keys()), sharex=True, figsize=(8, 10))

for press, data in output_dict.items():

    # write out data
    with open("trigonal_deformation_path_{}_GPa.txt".format(press), "w") as f:

        f.write("# Angle [°]; Optimum Volume [A**3]; Minimum Energy [meV]\n")

        for i in range(len(data[0])):

            f.write("{} {} {}\n".format(data[0][i], data[1][i], data[2][i]))

    # plotting
    # define index such that pressure increases from bottom to top
    index = -1-press_list.index(press)
    axes[index].plot(data[0], np.array(data[2])*1000-np.min(data[2])*1000, color="#DC143C", linewidth=2)
    axes[index].set_ylim(0, 100)
    axes[index].set_yticks(np.arange(0, 101, 25))
    axes[index].set_ylabel("{} GPa".format(press), fontsize=15)
    if press_list.index(press) == len(output_dict.keys())-1:  # Title for the top subplot
        axes[index].set_title("Relative Enthalpy [meV]", fontsize=15)
    axes[index].tick_params(axis='both', which='major', labelsize=12)

axes[-1].set_xticks(np.arange(60,111,5))
axes[-1].set_xlabel("Angle [°]", fontsize=15)
plt.savefig("trigonal_deformation_path_{}-{}_GPa.png".format(min_press, max_press), dpi=300)










