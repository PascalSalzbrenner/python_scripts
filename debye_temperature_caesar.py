# script to calculate the high-temperature-entropy Debye temperature from Caesar freq_dos.agr output
# Caesar uses atomic units so the frequencies are in Hartree
# implements equation (12) from Chen, Sundman - Acta Materialia Volume 49, Issue 6, 2 April 2001, Pages 947-961 to calculate the Debye frequency
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np

# set up conversion factor - 1/kB in K/Ha
conversion_factor = 315775.21536

# open dos input file
dosfile = open("freq-dos.agr", "r")

# the dos must be normalised to 1, so we divide the log moment by the integrated dos
# set up integrals
log_moment = 0
integrated_dos = 0

# the integral is ln(freq)*dos*d_freq
# we use the middle between two successive points
# have the freq and dos be initially undefined in order to check whether previous data exists
prev_freq = None
prev_dos = None

# read over file
for line in dosfile:

	data = line.split()
	new_freq = float(data[0])
	new_dos = float(data[1])

	if prev_freq:
		# only do this if we already have data from a previous step
		d_freq = new_freq - prev_freq
		freq = (new_freq+prev_freq)/2
		dos = (new_dos+prev_dos)/2

		log_moment += np.log(freq)*dos*d_freq
		integrated_dos += dos*d_freq

	prev_freq = new_freq
	prev_dos = new_dos

dosfile.close()

# calculate Debye frequency and temperature
debye_frequency = np.exp(1/3+log_moment/integrated_dos)
debye_temperature = conversion_factor * debye_frequency

# write output
with open("debye_temperature.dat", "w") as outfile:
	outfile.write("The Debye frequency of {} is {} Hartree.\n".format(seed, debye_frequency))
	outfile.write("The Debye temperature, therefore, is {} K.".format(debye_temperature))


