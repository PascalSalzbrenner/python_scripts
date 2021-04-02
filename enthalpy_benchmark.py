# script to calculate the difference between DFT and QMC calculations for several structures for enthalpy as a function of pressure
# currently rather specific, but can easily be generalised

# requires as input a QMC_enthalpy.dat file giving, well, the QMC enthalpies, in the format

# header
#
# # structure name
# pressure enthalpy
# ...

# and a DFT_enthalpy.dat doing the same for the DFT enthalpies
# denser pressure spacing in DFT is fine; the script will automatically find the DFT data point closest in pressure to the QMC data point

# calculates the mean error, equation 1 in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.184106
# and the root mean squared error, where the difference is squared and the square root of the sum is taken before dividing by the num_points

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

# energy units are assumed meV/atom, but this can easily be made an input variable
units = "meV/atom"

# set up data containers
# will contain the data in the form "structure_name": [[pressures], [enthalpies]]
qmc_dict = {}
dft_dict = {}

# will contain the data in the form "structure_name": [[differences], [squared_differences]]
diff_dict = {}

# will contain the data in the form "structure_name": [mean, rms]
means_dict = {}

# read QMC enthalpies
qmc_file = open("QMC_enthalpy.dat", "r")

# read past header
while qmc_file.readline().lstrip().startswith("#"):
    continue

# when this is done, we have reached the empty line indicating the end of the header - first structure name in the next line
for line in qmc_file:
    if len(line.split()) == 0:
        # empty line
        continue
    elif line.lstrip().startswith("#"):
        # structure name
        structure_name = line.lstrip().lstrip("#").lstrip().rstrip("\n")

        qmc_dict[structure_name] = [[],[]]
    else:
        # data point
        data = line.split()
        qmc_dict[structure_name][0].append(float(data[0]))
        qmc_dict[structure_name][1].append(float(data[1]))

qmc_file.close()

# initialise qmc_press_counter
qmc_press_counter = 1

# read DFT enthalpies
dft_file = open("DFT_enthalpy.dat", "r")

# read past header
while dft_file.readline().lstrip().startswith("#"):
    continue

# set flag to tell the program to read data
read = True

# when this is done, we have reached the empty line indicating the end of the header - first structure name in the next line
for line in dft_file:
    if len(line.split()) == 0:
        # empty line
        continue
    elif line.lstrip().startswith("#"):

        if qmc_press_counter > 1:
            # this is not the first structure, hence there are still variables from the last structure defined

            # check if the difference between the last DFT pressure read and the QMC pressure currently in the buffer is sufficiently small
            if prev_press_diff < 1:
                # the 1 is rather arbitrary, but seems sensible
                dft_dict[structure_name][0].append(prev_press)
                dft_dict[structure_name][1].append(prev_enthalpy)

        # structure name
        structure_name = line.lstrip().lstrip("#").lstrip().rstrip("\n")

        dft_dict[structure_name] = [[],[]]

        # from now on, read data again
        read = True

        # take the first QMC pressure for this structure
        qmc_press = qmc_dict[structure_name][0][0]
        # set prev_press_diff to a string to indicate we are dealing with a new structure
        prev_press_diff="press_diff"
        # set the qmc_press_counter to 1 - the index of the next pressure we have to read out
        qmc_press_counter=1
    else:
        # data point

        if read:
            data = line.split()

            press_diff=abs(float(data[0]) - qmc_press)

            # search for the DFT pressure closest to the current QMC pressure
            if prev_press_diff != "press_diff":
                # not a new structure, prev_press_diff will be a float
                if press_diff <= prev_press_diff:
                    # the new DFT pressure is just as close or closer to the QMC pressure as the previous one - continue searching
                    prev_press_diff=press_diff
                    prev_press=float(data[0])
                    prev_enthalpy=float(data[1])
                else:
                    # previous pressure was closer to the QMC pressure and hence is the best possible estimate
                    dft_dict[structure_name][0].append(prev_press)
                    dft_dict[structure_name][1].append(prev_enthalpy)

                    # read new QMC pressure and iterate counter
                    qmc_press = qmc_dict[structure_name][0][qmc_press_counter]
                    # check if there are still QMC pressures to be read
                    if len(qmc_dict[structure_name][0]) > qmc_press_counter + 1:
                        qmc_press_counter += 1
                        prev_press_diff = abs(float(data[0]) - qmc_press)
                    else:
                        # if this is false, we have reached the end of the pressure list for this structure
                        # set a flag telling the program not to read DFT data until it reaches the next structure
                        read = False
            else:
                # the first pressure for the new structure

                if qmc_press < float(data[0]):
                    # QMC pressure is smaller - find the first that is >= the DFT pressure
                    qmc_pressure_list = qmc_dict[structure_name][0][:]
                    for pressure in qmc_pressure_list:
                        if qmc_pressure_list[i] < float(data[0]):
                            # not there yet - reduce the lists in qmc_dict[structure_name] by that element
                            qmc_dict[structure_name][0] = qmc_dict[structure_name][0][1:]
                            qmc_dict[structure_name][1] = qmc_dict[structure_name][1][1:]
                        else:
                            qmc_press = pressure
                            break

                prev_press_diff=abs(float(data[0]) - qmc_press)
                prev_press=float(data[0])
                prev_enthalpy=float(data[1])

dft_file.close()

# the very last DFT value yet will not have been added to the dictionary - check if it should be
if prev_press_diff < 1:
    # for all the other structures, this is handled within the loop when the next structure begins
    dft_dict[structure_name][0].append(prev_press)
    dft_dict[structure_name][1].append(prev_enthalpy)

# calculate differences and averages
for key, value in dft_dict.items():

    # within the loop for reading the DFT data, it is ensured that there is never going to be less QMC than DFT data
    # here, we ensure the opposite
    if len(value[0]) < len(qmc_dict[key][0]):
        qmc_dict[key][0] = qmc_dict[key][0][:len(value[0])]
        qmc_dict[key][1] = qmc_dict[key][1][:len(value[1])]

    if key in qmc_dict.keys():
        diff_dict[key] = [np.array(np.array(value[1])-np.array(qmc_dict[key][1]))]
        diff_dict[key].append(np.power(diff_dict[key][0], 2))

        means_dict[key] = [np.mean(diff_dict[key][0]), np.sqrt(np.mean(diff_dict[key][1]))]

# write out data
outfile = open("enthalpy_differences.dat", "w")

for key, value in diff_dict.items():

    outfile.write("# {}\n".format(key))
    outfile.write("# QMC Pressure [GPa]; DFT Pressure [GPa]; H_DFT-H_QMC [{0}]; (H_DFT-H_QMC)**2 [{0}]\n".format(units))

    for i in range(len(value[0])):
        outfile.write("{} {} {} {}\n".format(qmc_dict[key][0][i], dft_dict[key][0][i], value[0][i], value[1][i]))

    outfile.write("\n")

outfile.write("# Structure; Mean [{0}]; RMS [{0}]\n".format(units))

# set up list for means and rms
means = []
rms = []

for key, value in means_dict.items():

    outfile.write("{} {} {}\n".format(key, value[0], value[1]))
    means.append(value[0])
    rms.append(value[1])

outfile.write("\n")

all_structures_mean = np.mean(np.array(means))
all_structures_rms = np.mean(np.array(rms))

outfile.write("# Mean difference between DFT and QMC, averaged over all structures: {} {}\n".format(all_structures_mean, units))
outfile.write("# RMS difference between DFT and QMC, averaged over all structures: {} {}\n".format(all_structures_rms, units))

outfile.close()
