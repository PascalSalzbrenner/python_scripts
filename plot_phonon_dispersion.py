# python postprocessing script to plot the phonon dispersions calculated by the in-house (bm418) Caesar code
# written by Pascal Salzbrenner, pts28

import re
import sys
import numpy as np
from math import ceil

class UnitError(Exception):
    """Exception raised when a conversion to something other than per_cm or meV is requested"""

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

# convert the frequencies given by Caesar in Hartree to cm^(-1) or meV

convert_to = sys.argv[1]

if convert_to.lower() == "per_cm":
    conversion_factor = 219474.63
elif convert_to.lower() == "mev":
    conversion_factor = 27211.399
else:
    raise UnitError("unit conversion",
    "You requested a conversion to {}, but the only possible options are per_cm or meV. Please enter one of the two.".format(convert_to))

dispersion_hartree = open("phonon_dispersion_curve.dat", "r")
dispersion_converted = open("phonon_dispersion_curve_{}.dat".format(convert_to), "w")

max_frequency = 0 # it will definitely be larger than this, as 0 is the minimum

for line in dispersion_hartree:
    frequencies = line.split()

    out_line = frequencies[0] # the value in k-space

    max_point_frequency = conversion_factor*np.double(frequencies[-1])

    if max_point_frequency > max_frequency:
        max_frequency = max_point_frequency

    for frequency in frequencies[1:]:
        out_line += "  {}".format(str(conversion_factor*np.double(frequency)))

    dispersion_converted.write("{}\n".format(out_line))

dispersion_hartree.close()
dispersion_converted.close()

ymax = round_up(max_frequency, 1)

# write gnuplot script

high_symmetry_points = open("high_symmetry_points.dat", "r")
high_symmetry_names = open("path.dat", "r")

plotfile = open("phonon_dispersion.gnu", "w")

hsp_name_list = []
hsp_location_list = []

# read high symmetry point names and locations

for hsp in high_symmetry_points:

    hsp_name = high_symmetry_names.readline().split()[-1].lstrip("#")
    hsp_location = hsp.split()[1]
    hsp_name_list.append(hsp_name)
    hsp_location_list.append(hsp_location)

high_symmetry_names.close()
high_symmetry_points.close()

plotfile.write("set terminal postscript eps colour\n")
plotfile.write("set style data points\n")
plotfile.write("set output '| epstopdf --filter --outfile=phonon_dispersion.pdf'\n")
plotfile.write("set mytics 2\n")

plotfile.write("set yrange [0:{}]\n".format(ymax))
plotfile.write("set xrange [0:{}]\n".format(hsp_location)) # the last hsp_location read is the highest value for which the phonon dispersions are calculated

xtics_string = "set xtics ("

for i in range(len(hsp_name_list)):

    # add lines at the hsp locations

    plotfile.write("set arrow from {0}, 0 to {0}, {1} nohead\n".format(hsp_location_list[i], ymax))

    # build xtics string
    if "G" in hsp_name_list[i].upper():
        xtics_string += "'{}' {}, ".format(re.sub(r"[Gg][(AMMA)(amma)]*", "{/Symbol G}", hsp_name_list[i]),  hsp_location_list[i])
    else:
        xtics_string += "'{}' {}, ".format(hsp_name_list[i], hsp_location_list[i])

xtics_string = xtics_string.rstrip(", ")

plotfile.write(xtics_string)

plotfile.write(")\n")

if convert_to.lower() == "per_cm":
    plotfile.write("set ylabel 'Wavenumber ~{/Symbol n}{.8\_} [cm^{-1}]'\n")
elif convert_to.lower() == "mev":
    plotfile.write("set ylabel 'Energy [meV]'\n")

plot_string = "plot "

for i in range(1, len(frequencies)):

    plot_string += "'phonon_dispersion_curve_{}.dat' u 1:{} pt 7 ps 0.3 lc rgb '#DC143C' notitle, ".format(convert_to, i+1)

plot_string = plot_string.rstrip(", ")
plot_string += "\n"

plotfile.write(plot_string)

plotfile.close()
