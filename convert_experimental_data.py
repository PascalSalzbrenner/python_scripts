# python postprocessing script to convert phonon dispersion data given in THz, cm^-1, or meV into each other
# phonon_experiment.dat's first non-comment line must give the units
# for cm^-1, both per_cm or cm^-1 are understood
# written by Pascal Salzbrenner, pts28

import numpy as np
from math import ceil

class UnitError(Exception):
    """Exception raised when a conversion to something other when the units of the input are not compatible"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

def round_up(num, num_int = None):
    """Function to round up to given integer place - eg 1 means tens (as 0 means single digits)
       If the decimal place is not given, round to the nearest integer whose only non-zero digit is the highest existing digit"""

    if not num_int:
        factor = 1
        while num > 10:
            num /= 10
            factor *= 10
        num = ceil(num)
        num *= factor
    else:
        num *= (10**(-num_int))
        num = ceil(num)
        num *= (10**num_int)

    return int(num)

convert_to = input("Which units would you like to convert to? (per_cm, meV, or THz) ")

infile = open("phonon_experiment.dat", "r")
outfile = open("phonon_experiment_{}.dat".format(convert_to), "w")

# read until first non-comment line, which gives the unit of the phonon data

for line in infile:

    if line.startswith("#"):
        # comment lines are written without processing
        outfile.write(line)
    else:
        convert_from = line.split()[0]
        outfile.write("{}\n".format(convert_to))
        break

# determine conversion factor

if convert_from.lower().startswith("t"):
    # THz
    if convert_to.lower().startswith("t"):
        conversion_factor = 1
    elif convert_to.lower().startswith("c") or convert_to.lower().startswith("p"):
        conversion_factor = 33.35640952
    elif convert_to.lower().startswith("m"):
        conversion_factor = 4.13567
    else:
        raise UnitError("Conversion from THz",
        "You requested conversion to {}, but this is not an option. Please enter either per_cm, cm^-1, meV, or THz".format(convert_to))
elif convert_from.lower().startswith("c") or convert_from.lower().startswith("p"):
    # cm^-1 or per_cm
    if convert_to.lower().startswith("t"):
        conversion_factor = 0.02998
    elif convert_to.lower().startswith("c") or convert_to.lower().startswith("p"):
        conversion_factor = 1
    elif convert_to.lower().startswith("m"):
        conversion_factor = 0.12399
    else:
        raise UnitError("Conversion from cm^-1",
        "You requested conversion to {}, but this is not an option. Please enter either per_cm, cm^-1, meV, or THz".format(convert_to))
elif convert_from.lower().startswith("m"):
    # meV
    if convert_to.lower().startswith("t"):
        conversion_factor = 0.24180
    elif convert_to.lower().startswith("c") or convert_to.lower().startswith("p"):
        conversion_factor = 8.06558
    elif convert_to.lower().startswith("m"):
        conversion_factor = 1
    else:
        raise UnitError("Conversion from meV",
        "You requested conversion to {}, but this is not an option. Please enter either per_cm, cm^-1, meV, or THz".format(convert_to))
else:
    raise UnitError("Conversion from {}".format(convert_from),
    "You requested conversion from {}, but this is not an option. Please ensure that the first non-comment line of your file is either per_cm, cm^-1, meV, or THz".format(convert_from))

for line in infile:
    frequencies = line.split()

    out_line = frequencies[0] # the value in k-space

    for frequency in frequencies[1:]:
        out_line += "  {}".format(str(conversion_factor*np.double(frequency)))

    outfile.write("{}\n".format(out_line))

infile.close()
outfile.close()
