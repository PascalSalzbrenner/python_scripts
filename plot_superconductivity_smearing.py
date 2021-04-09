# script to plot the output of the Quantum Espresso script lambda.x for different phonon grids
# any of lambda, omega_log, N(Ef), T_c can be plotted against degauss
# N will be the same for all phonon grid sizes as it only depends on the electron grid
# inputs: names of directories containing lambda_mu_<mu>.out files - must be given as command line arguments

# NB: I am not fully sure about the units as they're not given in the QE output, but the below are my best guesses

# written by Pascal Salzbrenner - pts28@cam.ac.uk

import sys
import math
import numpy as np

mu = input("At which mu* was the output calculated? ")
quantity = input("Which quantity would you like you plot (lambda, omega_log, N (eDOS at the Fermi level), Tc)? ")

if quantity.lower().startswith("l"):
    units = "dimensionless"
    column = "2"
elif quantity.lower().startswith("o"):
    units = "cm^{-1}"
    column = "3"
elif quantity.lower().startswith("n"):
    units = "electrons/cell"
    column = "4"
else:
    # Tc is default
    units = "K"
    column = "5"

plotfile = open("{}_degauss_mu_{}.gnu".format(quantity, mu), "w") # gnuplot file

plotfile.write("set terminal postscript eps enhanced colour\n")
plotfile.write("set style data linespoints\n")
plotfile.write("set output '| epstopdf --filter --outfile={}_degauss_mu_{}.pdf'\n".format(quantity, mu))

plotfile.write("set key top right\n")
plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")

# redefine gnuplot linetypes with nice colours
plotfile.write("set linetype 1 lc rgb '#DC143C' lt 1\n")
plotfile.write("set linetype 2 lc rgb '#D95F02' lt 1\n")
plotfile.write("set linetype 3 lc rgb '#E6AB02' lt 1\n")
plotfile.write("set linetype 4 lc rgb '#66A61E' lt 1\n")
plotfile.write("set linetype 5 lc rgb '#8000C4' lt 1\n")
plotfile.write("set linetype 6 lc rgb '#7570B3' lt 1\n")
plotfile.write("set linetype 7 lc rgb '#E7298A' lt 1\n")
plotfile.write("set linetype 8 lc rgb '#1E90FF' lt 1\n")
plotfile.write("set linetype 9 lc rgb '#1B9E77' lt 1\n")
plotfile.write("set linetype 10 lc rgb '#B8860B' lt 1\n")
plotfile.write("set linetype 11 lc rgb '#20C2C2' lt 1\n")
plotfile.write("set linetype cycle 11\n")

plotfile.write("set xlabel 'Smearing width [Ry]'\n")
plotfile.write("set ylabel '{} [{}]'\n".format(quantity, units))
plotfile.write("set mxtics 2\n")
plotfile.write("set mytics 2\n")

plot_string = "plot" # there is one plot command at the end, so we generate the string and write it to energy_volume.gnu at the end

for directory in sys.argv[1:]:

    clean_dirname = directory.rstrip("/")

    # process data from lambda.x into more easily plotted output
    lambda_file = open("{}/lambda_mu_{}.out".format(clean_dirname, mu), "r")
    processed_file = open("{}/lambda_mu_{}_output_plot.dat".format(clean_dirname, mu), "w")

    # set up containers
    smearings = []
    lambdas = []
    omegas = []
    dos = []
    tcs = []

    processed_file.write("Gaussian smearing [Ry]; lambda [dimensionless]; omega_log [cm^-1]; eDOS at the Fermi level [electrons/cell]; Tc [Ry]\n")

    for line in lambda_file:
        if "degauss" in line:
            data = line.split()

            smearings.append(data[-1])
            dos.append(data[-4])
        elif "lambda" in line:
            # this line divides the first el_ph_nsigma lines containing the degauss from the second set summarising lambda, omega_log, Tc
            # does not contain any data itself
            continue
        else:
            # line contains lambda, omega_log, Tc
            data = line.split()
            lambdas.append(data[0])
            omegas.append(data[1])
            tcs.append(data[2])

    lambda_file.close()

    for i in range(len(smearings)):
        processed_file.write("{} {} {} {} {}\n".format(smearings[i], lambdas[i], omegas[i], dos[i], tcs[i]))

    processed_file.close()

    plot_string += " '{0}/lambda_mu_{1}_output_plot.dat' u 1:{2} w linespoints pt 7 ps 1 title '{0}',".format(clean_dirname, mu, column)

plot_string=plot_string.rstrip(",")
plotfile.write(plot_string)

plotfile.close()
