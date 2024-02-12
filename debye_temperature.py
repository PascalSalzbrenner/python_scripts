# script to calculate the Debye temperature from wobble seed-dos.agr output
# implements equations 6.33 and 6.34 from G. Grimvall - Thermophysical Properties of Materials
# allows the calculation of the Debye frequency from arbitrary frequency moments
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np

# define function to integrate for the calculation of different moments
def moment_function(freq, n):
	"""The moment of the frequency which is used"""

	if n < -2:
		raise(ValueError, "Moments smaller than -2 cannot be calculated.")
	elif n == 0:
		return np.log(freq)
	else:
		return freq**n

# read seed
seed = sys.argv[1]

# read moment index
moment_index = int(sys.argv[2])

# open dos input file
dosfile = open("{}-dos.agr".format(seed), "r")

# the dos must be normalised to 1, so we divide the moment by the integrated dos
# set up integrals
frequency_moment = 0
integrated_dos = 0

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
			unit = data[4].lstrip('(').rstrip(')"')

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
			dos = (new_dos+prev_dos)/2

			frequency_moment += moment_function(freq, moment_index)*dos*d_freq
			integrated_dos += dos*d_freq

		prev_freq = new_freq
		prev_dos = new_dos

dosfile.close()

# calculate Debye frequency and temperature
if moment_index == 0:
	debye_frequency = np.exp(1/3+frequency_moment/integrated_dos)
else:
	debye_frequency = ((frequency_moment/integrated_dos)*(moment_index+3)/3)**(1/moment_index)

debye_temperature = conversion_factor * debye_frequency

# write output
with open("debye_temperature_{}_moment.dat".format(moment_index), "w") as outfile:
	outfile.write("The Debye frequency of {} is {} {}.\n".format(seed, debye_frequency, unit))
	outfile.write("The Debye temperature, therefore, is {} K.".format(debye_temperature))


