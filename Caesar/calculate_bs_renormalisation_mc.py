# script to calculate the average renormalisation of the difference between two bands due to electron-phonon coupling
# takes as input the raw data generated by the bs_mc.sh script of the in-house Caesar code

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np

############################################ define function to propogate errors when averaging ############################################

def propagate_error(uncertainties):
    """Function to propagate the uncertainties associated with quantities whose arithmetic mean is taken to that mean
    NB: assumes the quantities being averaged are uncorrelated, that is, their covariance is 0
    :param list uncertainties: list of floats expressing the uncertainties of the quantities being average

    :returns float uncertainty: float indicating the propagated uncertainty"""

    uncertainties = np.array(uncertainties)

    return (np.sqrt(np.sum(uncertainties**2)))/len(uncertainties)

############################################################# reading of input #############################################################

k_point = int(sys.argv[1]) # the k-point number (first number of the two in the output files), ordered as in IBZKPT
lower_band_index = int(sys.argv[2]) # the index of the lower band of the two bands whose difference we are calculating (eg the VB)
lower_band_degeneracy = int(sys.argv[3]) # the degeneracy of the band
higher_band_index = int(sys.argv[4]) # the index of the higher band of the two bands whose difference we are calculating (eg the CB)
higher_band_degeneracy = int(sys.argv[5]) # the degeneracy of the band
# indices refer to the band index in the corresponding static supercell calculation
# in the case of band degeneracy, they refer to the highest band

############################################# treatment of bs.k_point.band_index.dat data #############################################

# Average the band energies for the different degenerate bands, for both the lower and the higher band, and then calculate the difference
# lower band
for i in range(lower_band_degeneracy):
    bs_file = open("bs/bs.{}.{}.dat".format(k_point, lower_band_index-i), "r")
    lower_band_bs_temp = []

    for line in bs_file:
        lower_band_bs_temp.append(float(line))

    if i == 0:
        # store into the output vector
        lower_band_bs = np.array(lower_band_bs_temp)

    else:
        # store into temp vector for averaging with output vector
        lower_band_bs_temp = np.array(lower_band_bs_temp)
        lower_band_bs = (i*lower_band_bs + lower_band_bs_temp)/(i+1)

    bs_file.close()

# higher band
for i in range(higher_band_degeneracy):
    bs_file = open("bs/bs.{}.{}.dat".format(k_point, higher_band_index-i), "r")
    higher_band_bs_temp = []

    for line in bs_file:
        higher_band_bs_temp.append(float(line))

    if i == 0:
        # store into the output vector
        higher_band_bs = np.array(higher_band_bs_temp)
    else:
        # store into temp vector for averaging with output vector
        higher_band_bs_temp = np.array(higher_band_bs_temp)
        higher_band_bs = (i*higher_band_bs + higher_band_bs_temp)/(i+1)

    bs_file.close()

# generate output for both the averaged lower and higher sets of bands, as well as their difference
bs_lower_lines = []
bs_higher_lines = []
bs_diff_lines = []

# comment lines
bs_lower_lines.append("# Averaged value over all degenerate bands of band {} at k-point {}".format(lower_band_index, k_point))
bs_higher_lines.append("# Averaged value over all degenerate bands of band {} at k-point {}".format(higher_band_index, k_point))
bs_diff_lines.append("# Energy difference between bands {} and {} at k-point {}".format(lower_band_index,
                                  higher_band_index, k_point))

for i in range(len(lower_band_bs)):
    bs_lower_lines.append(str(lower_band_bs[i]))
    bs_higher_lines.append(str(higher_band_bs[i]))
    bs_diff_lines.append(str(higher_band_bs[i]-lower_band_bs[i]))

# write output files
with open("average_bs_{}_{}.dat".format(k_point, lower_band_index), "w") as bs_lower_outfile:
    bs_lower_outfile.write("\n".join(bs_lower_lines))

with open("average_bs_{}_{}.dat".format(k_point, higher_band_index), "w") as bs_higher_outfile:
    bs_higher_outfile.write("\n".join(bs_higher_lines))

with open("bs_difference_{}_{}_{}.dat".format(k_point, lower_band_index, higher_band_index), "w") as bs_average_outfile:
    bs_average_outfile.write("\n".join(bs_diff_lines))

############################################# treatment of average.k_point.band_index.dat data #############################################

# calculate running averages for the different degenerate bands, for both the lower and the higher band, and then calculate the difference
# lower band
for i in range(lower_band_degeneracy):
    running_average_file = open("bs/average.{}.{}.dat".format(k_point, lower_band_index-i), "r")
    lower_band_running_average_temp = []

    for line in running_average_file:
        lower_band_running_average_temp.append(float(line))

    if i == 0:
        # store into the output vector
        lower_band_running_average = np.array(lower_band_running_average_temp)

        # make container for all lower bands which are not degenerate at the point (even though they will be degenerate on average)
        non_degenerate_lower_bands = [lower_band_index]

        # store the running average in a separate vector to compare to the running average for the next band to determine degeneracy
        prev_running_average = np.array(lower_band_running_average_temp)
    else:
        # store into temp vector for averaging with output vector
        lower_band_running_average_temp = np.array(lower_band_running_average_temp)
        lower_band_running_average = (i*lower_band_running_average + lower_band_running_average_temp)/(i+1)

        # compare the running average for this band to that for the previous band
        if not np.isclose(prev_running_average, lower_band_running_average_temp).all():
            # non-degenerate
            non_degenerate_lower_bands.append(lower_band_index-i)
            # update running average to that of this band
            prev_running_average = np.copy(lower_band_running_average_temp)

    running_average_file.close()

# higher band
for i in range(higher_band_degeneracy):
    running_average_file = open("bs/average.{}.{}.dat".format(k_point, higher_band_index-i), "r")
    higher_band_running_average_temp = []

    for line in running_average_file:
        higher_band_running_average_temp.append(float(line))

    if i == 0:
        # store into the output vector
        higher_band_running_average = np.array(higher_band_running_average_temp)

        # make container for all higher bands which are not degenerate at the point (even though they will be degenerate on average)
        non_degenerate_higher_bands = [higher_band_index]

        # store the running average in a separate vector to compare to the running average for the next band to determine degeneracy
        prev_running_average = np.array(higher_band_running_average_temp)
    else:
        # store into temp vector for averaging with output vector
        higher_band_running_average_temp = np.array(higher_band_running_average_temp)
        higher_band_running_average = (i*higher_band_running_average + higher_band_running_average_temp)/(i+1)

        # compare the running average for this band to that for the previous band
        if not np.isclose(prev_running_average, higher_band_running_average_temp).all():
            # non-degenerate
            non_degenerate_higher_bands.append(higher_band_index-i)
            # update running average to that of this band
            prev_running_average = np.copy(higher_band_running_average_temp)

    running_average_file.close()

# generate output for both the averaged lower and higher sets of bands, as well as their difference
running_average_lower_lines = []
running_average_higher_lines = []
running_average_diff_lines = []

# comment lines
running_average_lower_lines.append("# Running average over all degenerate bands for band {} at k-point {}".format(lower_band_index, k_point))
running_average_higher_lines.append("# Running average over all degenerate bands for band {} at k-point {}".format(higher_band_index, k_point))
running_average_diff_lines.append("# Running average for the energy difference between bands {} and {} at k-point {}".format(lower_band_index,
                                  higher_band_index, k_point))

for i in range(len(lower_band_running_average)):
    running_average_lower_lines.append("{} {}".format(i+1, lower_band_running_average[i]))
    running_average_higher_lines.append("{} {}".format(i+1, higher_band_running_average[i]))
    running_average_diff_lines.append("{} {}".format(i+1, higher_band_running_average[i]-lower_band_running_average[i]))

# write output files
with open("average_band_{}_{}.dat".format(k_point, lower_band_index), "w") as running_average_lower_outfile:
    running_average_lower_outfile.write("\n".join(running_average_lower_lines))

with open("average_band_{}_{}.dat".format(k_point, higher_band_index), "w") as running_average_higher_outfile:
    running_average_higher_outfile.write("\n".join(running_average_higher_lines))

with open("average_band_difference_{}_{}_{}.dat".format(k_point, lower_band_index, higher_band_index), "w") as running_average_outfile:
    running_average_outfile.write("\n".join(running_average_diff_lines))

# write lines to plot running average
running_average_plotlines = ["set terminal postscript eps colour", "set style data linespoints",
                             "set output '| epstopdf --filter --outfile=running_average_{}_{}_{}.pdf'".format(k_point, lower_band_index,
                             higher_band_index), "unset xlabel", "set ylabel '{/Symbol D}E_{G} [eV]'",
                             "plot 'average_band_difference_{}_{}_{}.dat' u 1:2 w linespoints pt 7 ps 1.5 lc rgb '#DC143C' notitle".format(
                             k_point, lower_band_index, higher_band_index)]

# write file to plot running average
with open("running_average_{}_{}_{}.gnu".format(k_point, lower_band_index, higher_band_index), "w") as running_average_plotfile:
    running_average_plotfile.write("\n".join(running_average_plotlines))

################################ calculate the difference between the two bands at the static lattice level ################################

# open static EIGENVAL file
eigenval = open("configurations/static/EIGENVAL", "r")

# read past the 6 header lines
for line in eigenval:
    if not line.split():
        break

# read past k_point-1 blocks of eigenvalues to get to the block corresponding to k_point
counter = 1

while counter < k_point:
    for line in eigenval:
        if not line.split():
            counter += 1
            break

# read past the k-point coordinates
eigenval.readline()

# read to find the lower and higher band indices
for line in eigenval:
    data = line.split()
    if not data:
        # if an empty line is encountered, we can assume this means one of the two band indices wasn't found
        raise ValueError("One of the band indices you supplied does not refer to a band in this calculation.")

    if int(data[0]) == lower_band_index:
        lower_band_energy = float(data[1])
    elif int(data[0]) == higher_band_index:
        higher_band_energy = float(data[1])
        break
    else:
        continue

eigenval.close()

energy_difference = higher_band_energy - lower_band_energy

############################################## treatment of mean.k_point.band_index.dat data ##############################################

# set up arrays to contain the means and uncertainties for the different degenerate bands
lower_band_means = []
lower_band_uncertainties = []
higher_band_means = []
higher_band_uncertainties = []

# read the means and uncertainties of all different bands for both higher and lower sets
for band in non_degenerate_lower_bands:
    with open("bs/mean.{}.{}.dat".format(k_point, band), "r") as band_mean_file:
        mean, uncertainty = list(map(float, band_mean_file.readline().split()))
        lower_band_means.append(mean)
        lower_band_uncertainties.append(uncertainty)

for band in non_degenerate_higher_bands:
    with open("bs/mean.{}.{}.dat".format(k_point, band), "r") as band_mean_file:
        mean, uncertainty = list(map(float, band_mean_file.readline().split()))
        higher_band_means.append(mean)
        higher_band_uncertainties.append(uncertainty)

# calculate means and difference
lower_band_mean = np.mean(lower_band_means)
higher_band_mean = np.mean(higher_band_means)
mean_difference = higher_band_mean - lower_band_mean
renormalisation = mean_difference - energy_difference # change in the difference between the two bands caused by the phonons
# positive means electron-phonon coupling increases the energy difference

# calculate uncertainties
lower_band_uncertainty = propagate_error(lower_band_uncertainties)
higher_band_uncertainty = propagate_error(higher_band_uncertainties)
# the uncertainty of the difference between higher and lower band is also the uncertainty of the renormalisation
difference_uncertainty = np.sqrt(lower_band_uncertainty**2+higher_band_uncertainty**2)

# write output
with open("bs/degeneracy_average_{}_{}.dat".format(k_point, lower_band_index), "w") as outfile:
    outfile.write("# Average energy [eV]; Uncertainty [eV]\n")
    outfile.write("{} {}\n".format(lower_band_mean, lower_band_uncertainty))

with open("bs/degeneracy_average_{}_{}.dat".format(k_point, higher_band_index), "w") as outfile:
    outfile.write("# Average energy [eV]; Uncertainty [eV]\n")
    outfile.write("{} {}\n".format(higher_band_mean, higher_band_uncertainty))

with open("renormalisation_{}_{}_{}.dat".format(k_point, lower_band_index, higher_band_index), "w") as outfile:
    outfile.write("# Average energy [eV]; Uncertainty [eV]\n")
    outfile.write("# First line contains the difference between the bands, the second line the change compared to the static lattice\n")
    outfile.write("{} {}\n".format(mean_difference, difference_uncertainty))
    outfile.write("{} {}\n".format(renormalisation, difference_uncertainty))
