# script to plot the phonon densities of states from the following sources: wobble (EDDP), Caesar (DFT), as well as the associated Debye spectra
# Debye spectra are those associated with the high-temperature-entropy Debye frequency
# this is calculated from equation 6.34 of G. Grimvall - Thermophysical Properties of Materials to calculate the Debye frequency
# we normalise all the DOSes to 1
# output plot is always in meV

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np
import matplotlib.pyplot as plt

# define function to integrate DOS and log moment
def integrate_dos_log_moment(freq, dos):
    integrated_dos = 0
    log_moment = 0
    for i in range(len(freq)-1):
        d_freq = freq[i+1] - freq[i]
        freq_mean = (freq[i+1]+freq[i])/2
        dos_mean = (dos[i+1]+dos[i])/2

        integrated_dos += dos_mean*d_freq
        log_moment += np.log(freq_mean)*dos_mean*d_freq

    return integrated_dos, log_moment

# read seed
seed = sys.argv[1]

# set up log moment integrals for Debye frequency calculation
wobble_log_moment = 0
caesar_log_moment = 0

# set up dos integrals for normalisation
wobble_integrated_dos = 0
caesar_integrated_dos = 0

# set up data containers
wobble_freq = []
wobble_dos = []
caesar_freq = []
caesar_dos = []

##### wobble #####

# open dos input file
with open("{}-dos.agr".format(seed), "r") as dosfile:

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
                    # conversion factor is 1 since we are using meV
                    conversion_factor = 1
                elif unit == "THz":
                    conversion_factor = 0.24180
                elif unit == "cm-1":
                    conversion_factor = 8.06558

            else:
                continue

        elif "&" in line:
            # final line
            break
        else:
            data = line.split()
            wobble_freq.append(float(data[0])*conversion_factor)
            wobble_dos.append(float(data[1]))

# calculate and normalise DOS
wobble_freq = np.array(wobble_freq)
wobble_integrated_dos, wobble_log_moment = integrate_dos_log_moment(wobble_freq, wobble_dos)
wobble_dos = np.array(wobble_dos)/wobble_integrated_dos

# calculate Debye frequency
wobble_debye_frequency = np.exp(1/3+wobble_log_moment/wobble_integrated_dos)

##### Caesar #####

# open dos input file
with open("{}-freq_dos.dat".format(seed), "r") as dosfile:

    # Caesar gives energies in Hartree and we want to convert to meV
    conversion_factor = 27211.386

    # read over file
    for line in dosfile:

        data = line.split()
        caesar_freq.append(float(data[0])*conversion_factor)
        caesar_dos.append(float(data[1]))

        # only do this if we have at least two data points
        if len(caesar_freq) > 1:
            d_freq = caesar_freq[-1] - caesar_freq[-2]
            freq = (caesar_freq[-1]+caesar_freq[-2])/2
            dos = (caesar_dos[-1]+caesar_dos[-2])/2

            caesar_integrated_dos += dos*d_freq
            caesar_log_moment += np.log(freq)*dos*d_freq

# calculate and normalise DOS
caesar_freq = np.array(caesar_freq)
caesar_integrated_dos, caesar_log_moment = integrate_dos_log_moment(caesar_freq, caesar_dos)
caesar_dos = np.array(caesar_dos)/caesar_integrated_dos

# calculate Debye frequency
caesar_debye_frequency = np.exp(1/3+caesar_log_moment/caesar_integrated_dos)

##### plot #####

# define function for Debye spectra
def debye_dos(w, omega_D):
    if w <= omega_D:
        return 3*w**2/omega_D**3
    else:
        return 0

wobble_debye_dos = [debye_dos(w, wobble_debye_frequency) for w in wobble_freq]
caesar_debye_dos = [debye_dos(w, caesar_debye_frequency) for w in caesar_freq]

plt.plot(wobble_freq, wobble_dos, label="EDDP - SOC", color="#E6AB02", linestyle="solid")
plt.plot(wobble_freq, wobble_debye_dos, label="EDDP - Debye", color="#E6AB02", linestyle="dashed")

plt.plot(caesar_freq, caesar_dos, label="DFT - SOC", color="#66A61E", linestyle="solid")
plt.plot(caesar_freq, caesar_debye_dos, label="DFT - Debye", color="#66A61E", linestyle="dashed")

plt.xlabel("Frequency (meV)")
plt.ylabel("DOS (normalised)")

# hide y-axis labels
plt.yticks([])

plt.legend()
plt.savefig("{}-phdos.png".format(seed), dpi=500)
