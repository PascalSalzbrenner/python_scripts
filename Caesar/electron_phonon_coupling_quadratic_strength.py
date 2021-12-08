# script to calculate the strength of the band gap renormalisation due to individual phonon modes and write a gnuplot script to plot a heat
# map of contributions

# input files
# ibz.dat file giving the q-points at which the electron-phonon coupling is evaluated
# bg_correction_kp.dat files for both conduction and valence band

# output files
# vb_renorm_heat_map.gnu: gnuplot script to plot the strength of the valence band renormalisation
# cb_renorm_heat_map.gnu: gnuplot script to plot the strength of the conduction band renormalisation
# bg_renorm_heat_map.gnu: gnuplot script to plot the strength of the band gap renormalisation
# vb_renorm_heat_map.dat: valence band renormalisation data in the format required by a gnuplot heat map
# cb_renorm_heat_map.dat: conduction band renormalisation data in the format required by a gnuplot heat map
# bg_renorm_heat_map.dat: band gap renormalisation data in the format required by a gnuplot heat map

# written by Pascal Salzbrenner, pts28

import numpy as np
import matplotlib.pyplot as plt

# define conversion factor for eV to cm^(-1) - taken from Inguscio and Fallani - Atomic Physics: Precise Measurements and Ultracold Matter
conversion_factor = 8065.54429

reciprocal_lattice_vector_1 = int(input("What is the first reciprocal lattice vector along which you would like to plot the heat map (1,2,3)? "))
reciprocal_lattice_vector_2 = int(input("What is the second reciprocal lattice vector along which you would like to plot the heat map (1,2,3)? "))
bz_slice = float(input("At what level of the third reciprocal lattice vector would you like to take the cross section of the BZ? "))
blocks = input("Along which blocks, characterised by different values of the q-point in direction {}, would you like to plot the mode-resolved coupling strengths? ('a1 a2 a3 ...') ".format(reciprocal_lattice_vector_1))
num_atoms = int(input("How many atoms are in your primitive unit cell? "))
num_modes = 3*num_atoms

blocks = [int(i) for i in blocks.split()]

# determine the reciprocal lattice_vector along which not to plot
for i in range(1,4):
    if i != reciprocal_lattice_vector_1 and i != reciprocal_lattice_vector_2:
        reciprocal_lattice_vector_no_plot = i
    else:
        continue

# read total and mode-resolved renormalisations per k-point for VB and CB
vb_kpoint_total_renorm = []
cb_kpoint_total_renorm = []
vb_kpoint_mode_renorm = []
cb_kpoint_mode_renorm = []

# also read mode frequencies from frequency.dat files
kpoint_mode_frequencies = []
# the frequency of mode j at k-point i will be accessible as kpoint_mode_frequencies[i][j-1]

# set up high minimum strength to iteratively determine the actual minimum strength
cb_min_strength = np.inf
vb_min_strength = np.inf
bg_min_strength = np.inf

vb_renorm_file = open("bs_valence_band/bg_correction_kp.dat", "r")

for line in vb_renorm_file:
    mode_renorm = []
    mode_frequencies = []
    if line.lstrip().startswith("k-point"):
        # the next num_modes lines have the mode-resolved renormalisation strengths in order of increasing frequency
        for i in range(num_modes):
            mode_data = vb_renorm_file.readline().split()
            mode_renorm.append(float(mode_data[2]))

            # check if this strength is smaller than the minimum strength hitherto found
            if np.abs(float(mode_data[2])) < vb_min_strength and np.abs(float(mode_data[2])) > 0:
                vb_min_strength = np.abs(float(mode_data[2]))

            # read from corresponding frequency file - the frequencies are of course the same for CB and VB, so must only be read once
            with open("bs_valence_band/frequency.{}.{}.dat".format(mode_data[0], mode_data[1]), "r") as freqfile:
                mode_frequencies.append(float(freqfile.readline())*conversion_factor) # frequency given in eV and converted to cm^(-1)

        vb_kpoint_mode_renorm.append(mode_renorm[:])
        kpoint_mode_frequencies.append(mode_frequencies[:])
        # next line contains the sum
        vb_kpoint_total_renorm.append(float(vb_renorm_file.readline().split()[1]))

vb_renorm_file.close()

cb_renorm_file = open("bs_conduction_band/bg_correction_kp.dat", "r")

for line in cb_renorm_file:
    mode_renorm = []
    mode_frequencies = []
    if line.lstrip().startswith("k-point"):
        # the next num_modes lines have the mode-resolved renormalisation strengths in order of increasing frequency
        for i in range(num_modes):
            mode_data = cb_renorm_file.readline().split()
            mode_renorm.append(float(mode_data[2]))

            # check if this strength is smaller than the minimum strength hitherto found
            if np.abs(float(mode_data[2])) < cb_min_strength and np.abs(float(mode_data[2])) > 0:
                cb_min_strength = np.abs(float(mode_data[2]))

        cb_kpoint_mode_renorm.append(mode_renorm[:])
        # next line contains the sum
        cb_kpoint_total_renorm.append(float(cb_renorm_file.readline().split()[1]))

cb_renorm_file.close()

# turn the lists into np.arrays for easy manipulation
vb_kpoint_total_renorm = np.array(vb_kpoint_total_renorm)
cb_kpoint_total_renorm = np.array(cb_kpoint_total_renorm)
vb_kpoint_mode_renorm = np.array(vb_kpoint_mode_renorm)
cb_kpoint_mode_renorm = np.array(cb_kpoint_mode_renorm)

# calculate the band gap renormalisation
bg_kpoint_total_renorm = cb_kpoint_total_renorm - vb_kpoint_total_renorm
bg_kpoint_mode_renorm = cb_kpoint_mode_renorm - vb_kpoint_mode_renorm

for i in bg_kpoint_mode_renorm:
    for j in i:
        if np.abs(j) < bg_min_strength and np.abs(j) > 0:
            bg_min_strength = np.abs(j)

# read in and store q-points and their number in tuples of the form (reciprocal_coordinate_1, reciprocal_coordinate_2, reciprocal_coordinate_3, number)
q_point_num = 0
q_points = []

with open("ibz.dat", "r") as q_point_file:

    for line in q_point_file:
        q_point = line.split()
        q_points.append((float(q_point[0]), float(q_point[1]), float(q_point[2]), q_point_num))
        #q_points.append((float(q_point[0]), -float(q_point[1]), -float(q_point[2]), q_point_num))
        #q_points.append((float(q_point[0]), float(q_point[2]), float(q_point[1]), q_point_num))
        # add points according to symmetry operations 4 & 17 of the zincblende reciprocal lattice in the Bilbao Crystallographic Server
        q_point_num += 1

# order the list of q_points in ascending order of the q_x coordinate
# for each q_x coordinate, ordered in ascending order of the q_y coordinate
# for each q_y coordinate, ordered in ascending order of the q_z coordinate
# this is the default of the internal sort algorithm
q_points.sort()

# write out the data files for use with the gnuplot scripts
# scale strengths by 1000
vb_data_file = open("vb_renorm_heat_map.dat", "w")

# determine the q-vectors corresponding to the blocks for mode-resolved plotting
counter = 1
mode_resolved_plot_q_points = []
plot_q_points = []

for i in range(len(q_points)):
    if np.isclose(q_points[i][reciprocal_lattice_vector_no_plot-1], bz_slice):
        vb_data_file.write("{} {} {}\n".format(q_points[i][reciprocal_lattice_vector_1-1], q_points[i][reciprocal_lattice_vector_2-1],
        1000*vb_kpoint_total_renorm[q_points[i][3]]))

        if counter in blocks:
            plot_q_points.append(q_points[i])

        if i == len(q_points)-1 or q_points[i+1][reciprocal_lattice_vector_1-1] != q_points[i][reciprocal_lattice_vector_1-1]:
            # a block is finished so Gnuplot requires a blank line
            vb_data_file.write("\n")

            if plot_q_points:
                mode_resolved_plot_q_points.append(plot_q_points[:])
                plot_q_points = []

            counter += 1

vb_data_file.close()

cb_data_file = open("cb_renorm_heat_map.dat", "w")

for i in range(len(q_points)):
    if np.isclose(q_points[i][reciprocal_lattice_vector_no_plot-1], bz_slice):
        cb_data_file.write("{} {} {}\n".format(q_points[i][reciprocal_lattice_vector_1-1], q_points[i][reciprocal_lattice_vector_2-1],
        1000*cb_kpoint_total_renorm[q_points[i][3]]))
        if i == len(q_points)-1 or q_points[i+1][reciprocal_lattice_vector_1-1] != q_points[i][reciprocal_lattice_vector_1-1]:
            # a block is finished so Gnuplot requires a blank line
            cb_data_file.write("\n")

cb_data_file.close()

bg_data_file = open("bg_renorm_heat_map.dat", "w")

for i in range(len(q_points)):
    if np.isclose(q_points[i][reciprocal_lattice_vector_no_plot-1], bz_slice):
        bg_data_file.write("{} {} {}\n".format(q_points[i][reciprocal_lattice_vector_1-1], q_points[i][reciprocal_lattice_vector_2-1],
        1000*bg_kpoint_total_renorm[q_points[i][3]]))
        if i == len(q_points)-1 or q_points[i+1][reciprocal_lattice_vector_1-1] != q_points[i][reciprocal_lattice_vector_1-1]:
            # a block is finished so Gnuplot requires a blank line
            bg_data_file.write("\n")

bg_data_file.close()

################################################# mode-resolved part along lines in the BZ #################################################

# set the x- and y-coordinates for the spaghetti plot, as well as strengths and colours and xticks with labels
x = []
frequencies = []
strengths_vb = []
strengths_cb = []
strengths_bg = []
colours_vb = []
colours_cb = []
colours_bg = []
xticks = []
xtick_labels = []

for i in range(len(mode_resolved_plot_q_points)):

    xticks.append(i)
    if i == 0:
        xtick_labels.append("[0.00 {:.2f} -0.50]".format(mode_resolved_plot_q_points[0][0][1]))
    else:
        xtick_labels.append("[0.00 {:.2f} 0.50] | [0.00 {:.2f} -0.50]".format(mode_resolved_plot_q_points[i-1][-1][1], mode_resolved_plot_q_points[i][0][1]))

    for j in range(len(mode_resolved_plot_q_points[i])):
        for k in range(num_modes):
            # for each point, there are num_modes modes, which all have the same x-value
            # reciprocal_lattice_vector_2 will always go from -0.5 to 0.5
            x.append(mode_resolved_plot_q_points[i][j][reciprocal_lattice_vector_2-1] + i+0.5)
            frequencies.append(kpoint_mode_frequencies[mode_resolved_plot_q_points[i][j][3]][k])

            # scale strengths by the minimum strength so that all strengths are larger than 1
            # then, take the logarithm for the plotting
            vb_scaled_strength=np.abs(vb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k])/(200*vb_min_strength)

            if vb_scaled_strength == 0:
                strengths_vb.append(10)
            else:
                strengths_vb.append(vb_scaled_strength)

            if vb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] < 0:
                colours_vb.append("b")
            elif vb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] > 0:
                colours_vb.append("r")
            else:
                colours_vb.append("black")

            cb_scaled_strength=np.abs(cb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k])/(200*cb_min_strength)

            if cb_scaled_strength == 0:
                strengths_cb.append(10)
            else:
                strengths_cb.append(cb_scaled_strength)

            if cb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] < 0:
                colours_cb.append("b")
            elif cb_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] > 0:
                colours_cb.append("r")
            else:
                colours_cb.append("black")

            bg_scaled_strength=np.abs(bg_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k])/(200*bg_min_strength)

            if bg_scaled_strength == 0:
                strengths_bg.append(10)
            else:
                strengths_bg.append(bg_scaled_strength)

            if bg_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] < 0:
                colours_bg.append("b")
            elif bg_kpoint_mode_renorm[mode_resolved_plot_q_points[i][j][3]][k] > 0:
                colours_bg.append("r")
            else:
                colours_bg.append("black")

xticks.append(len(mode_resolved_plot_q_points))
xtick_labels.append("[0.00 {:.2f} 0.50]".format(mode_resolved_plot_q_points[-1][-1][1]))

vb_fig = plt.figure()
vb_ax = vb_fig.add_subplot(111)
vb_ax.scatter(x, frequencies, s=strengths_vb, c=colours_vb)
vb_ax.set_xlabel("q-point [reciprocal coordinates]")
vb_ax.set_ylabel("Phonon wavenumber $[cm^{-1}]$")
vb_ax.set_xticks(xticks)
vb_ax.set_xticklabels(xtick_labels)
vb_fig.savefig("vb_dispersion_mode_strength.pdf")

cb_fig = plt.figure()
cb_ax = cb_fig.add_subplot(111)
cb_ax.scatter(x, frequencies, s=strengths_cb, c=colours_cb)
cb_ax.set_xlabel("q-point [reciprocal coordinates]")
cb_ax.set_ylabel("Phonon wavenumber $[cm^{-1}]$")
cb_ax.set_xticks(xticks)
cb_ax.set_xticklabels(xtick_labels)
cb_fig.savefig("cb_dispersion_mode_strength.pdf")

bg_fig = plt.figure()
bg_ax = bg_fig.add_subplot(111)
bg_ax.scatter(x, frequencies, s=strengths_bg, c=colours_bg)
bg_ax.set_xlabel("q-point [reciprocal coordinates]")
bg_ax.set_ylabel("Phonon wavenumber $[cm^{-1}]$")
bg_ax.set_xticks(xticks)
bg_ax.set_xticklabels(xtick_labels)
bg_fig.savefig("bg_dispersion_mode_strength.pdf")
