# script to calculate the high-temperature-entropy Debye temperature from wobble seed-dos.agr output
# implements equation (12) from Chen, Sundman - Acta Materialia Volume 49, Issue 6, 2 April 2001, Pages 947-961 to calculate the Debye frequency
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np

# read seed
seed = sys.argv[1]

# open dos input file
dosfile = open("{}-dos.agr".format(seed), "r")

# set up integral
log_moment = 0

# the integral is ln(freq)*dos*d_freq
# we use the middle between two successive points
# have the freq and dos be initially undefined in order to check whether previous data exists
prev_freq = None
prev_dos = None

# read over file
for line in dosfile:

	if line.startswith("@"):

		# these are comment lines - the only one which is relevant for us is the one telling us the frequency units
		if 'xaxis' in line and 'label' in line and '"' in line:
			# we parse for " as it is the only thing differentiating the line with the frequency units from the line with the font
			data = line.split()
			unit = line[4].lstrip('(').rstrip(')"')

			# determine conversation factor
			if unit == "meV":
				# conversion factor is 1/k_B in K/meV
				conversion_factor = 11.60451812
			elif unit == "THz":
				# conversion factor is h/k_B in K/THz
				conversion_factor = 47.99243
			elif unit == "cm-1":
				# conversion factor is hc/k_B in K cm
				conversion_factor = 1.4387769

		else:
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
			freq = (new_freq+prev_freq)/2
			dos = (new_dos+prev_freq)/2

			log_moment += np.log(freq)*dos*d_freq

		prev_freq = new_freq
		prev_dos = new_dos

dosfile.close()

# calculate Debye frequency and temperature
debye_frequency = np.exp((1+log_moment)/3)
debye_temperature = conversion_factor * debye_frequency

# write output
with open("debye_temperature.dat", "w") as outfile:
	outfile.write("The Debye temperature of {} is {} K.".format(seed, debye_temperature))
