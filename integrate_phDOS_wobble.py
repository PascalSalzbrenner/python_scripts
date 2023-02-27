# script to integrate the phonon DOS from wobble seed-dos.agr output
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np

# read seed
seed = sys.argv[1]

# open dos input file
dosfile = open("{}-dos.agr".format(seed), "r")

# set up integral
integrated_dos = 0

# we use the middle between two successive points for the integral
# have the freq and dos be initially undefined in order to check whether previous data exists
prev_freq = None
prev_dos = None

# read over file
for line in dosfile:

	if line.startswith("@"):
			continue

	elif "&" in line:
		# final line
		break

	else:
		data = line.split()
		new_freq = float(data[0])
		new_dos = float(data[1])

		if prev_freq:
			# only do this if we already have data from a previous step
			d_freq = new_freq - prev_freq
			dos = (new_dos+prev_freq)/2

			integrated_dos += dos*d_freq

		prev_freq = new_freq
		prev_dos = new_dos

dosfile.close()

# write output
with open("integrated_dos.dat", "w") as outfile:
	outfile.write("The integral of the {} phonon DOS is {}.".format(seed, integrated_dos))