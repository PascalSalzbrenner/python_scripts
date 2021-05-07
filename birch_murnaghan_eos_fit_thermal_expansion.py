# script to fit to the Birch-Murnaghan equation of state (DOI: 10.1103/PhysRev.71.809) to the total free energy at a range of temperatures
# with the phonon contribution calculated by the groups in-house Caesar code
# inputs: energy.txt file as produced by my geometry_optimisation_cubic routine, interpolated_free_energy.dat file from caesar lte_harmonic routine
# run in thermal_expansion directory containing a static_lattice and a vibrating_lattice directory

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

# reading information about at which temperatures to evaluate phonon energies from command-line arguments

# NB: only works with integer steps in temperature. I think, in principle, Caesar allows non-integer steps in temperature, but I don't
# expect to need this (I don't think many people will at all), so it's not implemented here, as this makes my life a lot easier
lowest_T = int(sys.argv[1]) # the starting temperature
highest_T = int(sys.argv[2]) # the final temperature
T_step = int(sys.argv[3]) # the temperature step - must be N*dtemperature, where N is an integer larger than or equal to 1, and dtemperature
# is the temperature difference between two consecutive temperatures in interpolated_free_energy.dat

units = sys.argv[4] # can be eV or meV -  meV is the default, as it is the energy scale at which phonon energies play out

structure_type = sys.argv[5]

if len(sys.argv) == 7:
    plot_step = int(sys.argv[6]) # the temperature step between two plotted curves - can be left empty, in which case all temperatures at which the fit is calculated are plotted
    #must be an integer mutiple of T_step (no error if it isn't but nothing will be plotted)
else:
    plot_step = None

if units.lower() == "ev":
    factor = 1
else:
    factor = 1000

# determine prefactor according to structure_type
if structure_type.lower().startswith("p") or structure_type.lower().startswith("s"):
    # simple cubic
    prefactor = 1
elif structure_type.lower().startswith("f"):
    # FCC
    prefactor = 0.5
elif structure_type.lower().startswith("b"):
    # BCC
    prefactor = 0.25

# define x and y data
lattice_parameters = []
volumes = []
energies = {"static": []}

# read lattice constant, volume, and energy out of the energy.txt file
with open("static_lattice/energy.txt", "r") as energyfile:

    for line in energyfile:
        data = line.split()
        lattice_parameters.append(data[0]) # lattice constant
        volumes.append(float(data[1])) # volume
        energies["static"].append(float(data[2])*factor)

# create dictionary entries for all temperature the user has specified
for t in range(lowest_T, highest_T+1, T_step):
    energies["T_{}K".format(t)] = []

# read phonon energy information from interpolated_free_energy.dat
for lattice_parameter in lattice_parameters:

    with open("vibrating_lattice/lattice_parameter_{}/lte/interpolated_free_energy.dat".format(lattice_parameter), "r") as energyfile:

        for line in energyfile:
            data = line.split()

            if int(float(data[0])) >= lowest_T and int(float(data[0])) <= highest_T:
                if "T_{}K".format(int(float(data[0]))) in energies.keys():
                    energies["T_{}K".format(int(float(data[0])))].append(float(data[2])*factor
                    +energies["static"][lattice_parameters.index(lattice_parameter)])
                    # adding the static energy - lattice_parameters.index[lattice_parameter] is unique as each lattice_parameter is different
            elif int(float(data[0])) < lowest_T:
                continue
            elif int(float(data[0])) > highest_T:
                break

# write lattice parameters, volumes, and the corresponding energies
with open("energy_with_phonons.txt", "w") as outfile:

    outfile.write("#lattice parameter [A], Volume [A^3], static lattice,")

    T_string = ""

    for t in range(lowest_T, highest_T+1, T_step):
        T_string += " {}K,".format(t)

    T_string = T_string.rstrip(",")
    outfile.write("{}\n".format(T_string))

    for i in range(len(lattice_parameters)):

        outfile.write("{} {}".format(lattice_parameters[i], volumes[i]))

        energy_string = ""

        for energy in energies.values():
            energy_string += " {}".format(energy[i])
        outfile.write("{}\n".format(energy_string))

# renormalise the plots to the energy at the shortest lattice parameter at the highest temperature considered by subtracting this variable
renorm_factor = energies["T_{}K".format(highest_T)][0]

paramfile = open("birch_murnaghan_eos_parameters.txt", "w")
minfile = open("minimum_energy.txt", "w") # file containing the minimum energy at each temperature and the corresponding volume and lattice parameter
plotfile = open("energy_volume_temperature.gnu", "w") # gnuplot file to plot fit

minfile.write("#T [K], a [A], V_0 [A^3], E_0 [{}]\n".format(units))

plotfile.write("set terminal postscript eps colour font 'Helvetica,20'\n")
plotfile.write("set style data points\n")
plotfile.write("set output '| epstopdf --filter --outfile=energy_volume_temperature.pdf'\n")

if plot_step:
    plotfile.write("set key top right\n")
    plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")

# redefine gnuplot linetypes with nice colours
plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
plotfile.write("set linetype 2 lc rgb '#DC143C'\n")
plotfile.write("set linetype 3 lc rgb '#D95F02'\n")
plotfile.write("set linetype 4 lc rgb '#D95F02'\n")
plotfile.write("set linetype 5 lc rgb '#E6AB02'\n")
plotfile.write("set linetype 6 lc rgb '#E6AB02'\n")
plotfile.write("set linetype 7 lc rgb '#66A61E'\n")
plotfile.write("set linetype 8 lc rgb '#66A61E'\n")
plotfile.write("set linetype 9 lc rgb '#8000C4'\n")
plotfile.write("set linetype 10 lc rgb '#8000C4'\n")
plotfile.write("set linetype 11 lc rgb '#7570B3'\n")
plotfile.write("set linetype 12 lc rgb '#7570B3'\n")
plotfile.write("set linetype 13 lc rgb '#E7298A'\n")
plotfile.write("set linetype 14 lc rgb '#E7298A'\n")
plotfile.write("set linetype 15 lc rgb '#1E90FF'\n")
plotfile.write("set linetype 16 lc rgb '#1E90FF'\n")
plotfile.write("set linetype 17 lc rgb '#1B9E77'\n")
plotfile.write("set linetype 18 lc rgb '#1B9E77'\n")
plotfile.write("set linetype 19 lc rgb '#B8860B'\n")
plotfile.write("set linetype 20 lc rgb '#B8860B'\n")
plotfile.write("set linetype 21 lc rgb '#20C2C2'\n")
plotfile.write("set linetype 22 lc rgb '#20C2C2'\n")
plotfile.write("set linetype cycle 22\n")

plotfile.write("set xlabel 'Volume [A^3]'\n")
plotfile.write("set ylabel 'Energy [{}]'\n".format(units))
# by default energy is plotted against volume, but given that the lattice parameter is written to the energy_with_phonons file,
# this can easily be changed if desired
plotfile.write("set xrange [{}:{}]\n".format(math.floor(min(volumes)), math.ceil(max(volumes))))

plot_string = "plot" # there is one plot command at the end, so we generate the string and write it to energy_volume.gnu at the end

# fit Birch-Murnaghan EOS at every temperature
for temperature, energy in energies.items():

    parameters, covariance = curve_fit(birch_murnaghan_eos, volumes, energy, bounds=((2*min(energy), volumes[0], 0, -np.inf),
    (max(energy), volumes[-1], np.inf, np.inf)), maxfev=100000) # impose very general boundaries to ensure correct result
    # if the sampling in geometry_optimisation_cubic is anywhere near reasonable, this reliably gives the right result

    # write Birch-Murnaghan parameters at this temperature
    paramfile.write("Temperature: {}\n".format(temperature))
    paramfile.write("E_0 = {} {}\n".format(parameters[0], units))
    paramfile.write("V_0 = {} A^3\n".format(parameters[1]))
    paramfile.write("lattice parameter a = {} A\n".format((prefactor*np.abs(parameters[1]))**(1.0/3.0)))
    paramfile.write("B_0 = {} {}/A^3 = {} GPa\n".format(parameters[2], units, parameters[2]*160.21766208/factor)) #conversion (m)eV/A^3 to GPa
    paramfile.write("B_prime = {} [dimensionless]\n\n".format(parameters[3]))

    minfile.write("{} {} {} {}\n".format(temperature, (prefactor*np.abs(parameters[1]))**(1.0/3.0), parameters[1], parameters[0]))

    # write temperature-dependent variables
    if plot_step:
        if temperature != "static" and int(temperature.split("_")[1].rstrip("K"))%plot_step == 0:
            plotfile.write("E_0_{} = {}\n".format(temperature, parameters[0]-renorm_factor))
            plotfile.write("V_0_{} = {}\n".format(temperature, parameters[1]))
            plotfile.write("B_0_{} = {}\n".format(temperature, parameters[2]))
            plotfile.write("B_prime_{} = {}\n".format(temperature, parameters[3]))
            plotfile.write("birch_murnaghan_eos_{0}(x) = E_0_{0} + 9.0*V_0_{0}*B_0_{0}/(16.0)*(((V_0_{0}/x)**(2.0/3.0)-1)**(3)*B_prime_{0} + ((V_0_{0}/x)**(2.0/3.0)-1)**(2)*(6-4*(V_0_{0}/x)**(2.0/3.0)))\n".format(temperature))

            plot_string += " 'energy_with_phonons.txt' u 2:(${}-{}) w points pt 7 ps 1.5 title '{}', birch_murnaghan_eos_{}(x) lw 2 notitle,".format(list(energies.keys()).index(temperature)+3, renorm_factor, temperature.split("_")[1], temperature)
    else:
        plotfile.write("E_0_{} = {}\n".format(temperature, parameters[0]-renorm_factor))
        plotfile.write("V_0_{} = {}\n".format(temperature, parameters[1]))
        plotfile.write("B_0_{} = {}\n".format(temperature, parameters[2]))
        plotfile.write("B_prime_{} = {}\n".format(temperature, parameters[3]))
        plotfile.write("birch_murnaghan_eos_{0}(x) = E_0_{0} + 9.0*V_0_{0}*B_0_{0}/(16.0)*(((V_0_{0}/x)**(2.0/3.0)-1)**(3)*B_prime_{0} + ((V_0_{0}/x)**(2.0/3.0)-1)**(2)*(6-4*(V_0_{0}/x)**(2.0/3.0)))\n".format(temperature))

        plot_string += " 'energy_with_phonons.txt' u 2:(${}-{}) w points pt 7 ps 1.5 notitle, birch_murnaghan_eos_{}(x) lw 2 notitle,".format(list(energies.keys()).index(temperature)+3, renorm_factor, temperature)

if plot_step:
    plot_string += " 'minimum_energy.txt' u 3:($4-{}) every {}::1 w points pt 9 ps 1.5 lc rgb 'black' notitle\n".format(renorm_factor, plot_step/T_step)
else:
    plot_string += " 'minimum_energy.txt' u 3:($4-{}) w points pt 9 ps 1.5 lc rgb 'black' notitle\n".format(renorm_factor)
plotfile.write(plot_string)

paramfile.close()
minfile.close()
plotfile.close()
