# script to calculate the Debye temperature from Caesar freq_dos.dat output
# Caesar uses atomic units so the frequencies are in Hartree
# implements equations 6.33 and 6.34 from G. Grimvall - Thermophysical Properties of Materials
# allows the calculation of the Debye frequency from arbitrary frequency moments
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

# define function to integrate for the calculation of different moments
def moment_function(freq, n):
	"""The moment of the frequency which is used"""

	if n < -3:
		raise(ValueError, "Moments smaller than -3 cannot be calculated.")
	elif n == 0:
		return np.log(freq)
	else:
		return freq**n

# set up conversion factor - 1/kB in K/Ha
conversion_factor = 315775.21536

# open dos input file
dosfile = open("freq_dos.dat", "r")

# the dos must be normalised to 1, so we divide the log moment by the integrated dos
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
	outfile.write("The Debye frequency is {} Hartree.\n".format(debye_frequency))
	outfile.write("The Debye temperature, therefore, is {} K.".format(debye_temperature))


