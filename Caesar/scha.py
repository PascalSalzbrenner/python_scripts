# script to run an SCHA calculation to correct a harmonic mode in Caesar
# requires a mapping.dat file and an energy.dat file present in the directory
# mapping.dat file has the traditional Caesar format
# energy.dat has two columns where the first column is the displacement index (which must be converted into a displacement amplitude), and the second column the energy
# frequency in atomic units of the relevant mode must be give as input parameter

# the process has two parts: 
# first, find the parameters of the quartic fit to the actual energy landscape
# second, find the best self-consistent harmonic approximation to that quartic fit using the Newton-Raphson method
# self-consistent equations:
# f(omega) = omega**2 - omega2 - 3*lambda/(2*omega)*coth(omega/(2*kb*temp))
# f'(omega) = 2*omega + 3*lambda/(2*omega**2)*coth(omega/(2*kb*temp))+3*lambda/(4*omega*kb*temp)*cosh^2(omega/(2*kb*temp))

# based on the SCHA implemention of Bartomeu Monserrat
# references:
# SCHA idea: https://www.tandfonline.com/doi/abs/10.1080/14786440408520575
# equation for the self-consistent harmonic frequency: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.155411
# implementation: https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.1.033181

# script written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

######################################################## Setup ########################################################

def quartic_function(x, omega2, lam):
	"""The function fit to the data"""

	return 0.5*omega2*x*x+0.25*lam*x**4

def scha_newton_raphson(temp):
	"""Newton-Raphson solver for the self-consistent harmonic fit"""

	xn=np.sqrt(abs(omega2))
	for i in range(1000):
		omega=xn
		arg=omega/(2*kb*temp)
		f = omega**2 - omega2 - 3*lam/(2*omega)*(np.exp(arg)+np.exp(-arg))/(np.exp(arg)-np.exp(-arg))
		# arg=omega/(2*kb*temp)
		fd = 2*omega + 3*lam/(2*omega**2)*(np.exp(arg)+np.exp(-arg))/(np.exp(arg)-np.exp(-arg))+3*lam/(4*
			 omega*kb*temp)*(2/(np.exp(arg)-np.exp(-arg)))*(2/(np.exp(arg)-np.exp(-arg)))

		# set up next step
		xn1=xn-f/fd
    	
		if abs((xn1-xn)/xn) < tol:
			return xn1

		xn=xn1
    
	# did not converge in 1000 iterations
	print("Warning: Newton-Raphson self-consistent harmonic fit did not converge within 1000 iterations")
	return xn1

# constants
# divide by this number to convert from eV to Hartree (atomic units)
ev_to_hartree = 27.211399
# Boltzmann constant [Hartree/K]
kb=8.6173303E-5/ev_to_hartree
# Newton-Raphson tolerance
tol=1E-6

# input parameters
frequency = float(sys.argv[1])

# divide energy by the size of the supercell in which it was calculated
supercell_size = int(sys.argv[2])

################################################ carry out quartic fit ################################################

# determine displacement amplitude
# note that the amplitude is set by the first line in mapping.dat
# this amplitude is reached when the point's displacement index equals the maximum point set by the second line in mapping.dat
with open("mapping.dat", "r") as mapping_file:
	amp_index = float(mapping_file.readline().split()[0])
	max_point = float(mapping_file.readline().split()[1])

amplitude = amp_index / np.sqrt(2*np.abs(frequency)) # this is how it is defined in Caesar

# read energy file
disp_indices = []
energies = []

with open("energy.dat", "r") as energy_file:

	for line in energy_file:

		if line.startswith("#"):
			# skip comment lines
			continue

		disp_index, energy = line.split()
		disp_indices.append(float(disp_index))
		energies.append(float(energy))

		if np.isclose(float(disp_index), 0):
			static_energy = float(energy)

disp_amplitudes = np.array(disp_indices)*amplitude/max_point
energies = (np.array(energies)-static_energy)/(supercell_size*ev_to_hartree)

parameters, covariance = curve_fit(quartic_function, disp_amplitudes, energies)

omega2 = parameters[0]
lam = parameters[1]

displacement_fit_displacements = np.linspace(np.min(disp_amplitudes), np.max(disp_amplitudes), num=int((np.max(disp_indices)-np.min(disp_indices))*100))
displacement_fit_energies = quartic_function(displacement_fit_displacements, omega2, lam)

plt.plot(disp_amplitudes, energies, 'o', displacement_fit_displacements, displacement_fit_energies, '-')
plt.savefig("quartic_fit.png", dpi=300)
plt.close()

############################################ Self-consistent harmonic fit #############################################

temperatures = []
scha_energies_hartree = []

for temp in range(10, 411, 25):

	temperatures.append(temp)
	scha_energies_hartree.append(scha_newton_raphson(temp))

scha_energies_hartree = np.array(scha_energies_hartree)
scha_energies_ev = scha_energies_hartree * ev_to_hartree

with open("scha.dat", "w") as outfile:

	outfile.write("# Temperature [K]; Energy [Hartree]; Energy [eV]\n")

	for i in range(len(temperatures)):
		outfile.write("{} {} {}\n".format(temperatures[i], scha_energies_hartree[i], scha_energies_ev[i]))


















