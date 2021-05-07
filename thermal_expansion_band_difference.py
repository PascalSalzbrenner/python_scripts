# Python script which writes gnuplot files for plotting the (change in) band gap as a function of temperature

# input
# band_difference_temperature.dat file, with the format T [K], band gap [eV], change in band gap compared to static [meV]
# and with the correct k-point and band index in the name

# output
# band_gap_temperature.gnu, to plot the absolute value of the band gap as a function of temperature
# band_gap_change_temperature.gnu, to plot the change in the band gap compared to the static result as a function of temperature

import sys

lower_band_index = int(sys.argv[1]) # the index of the band whose energy we want to evaluate
k_point = sys.argv[2].split() # The k-point at which you want to evaluate the band energy, in the format kx ky kz

# create lists of lines for band gap gnuplot files
bd_plot_lines = ["set terminal postscript eps colour font 'Helvetica,20'", "set style data points",
"set output '| epstopdf --filter --outfile=bands_{}_{}_k_{}_{}_{}_difference_temperature.pdf'".format(lower_band_index, lower_band_index+1,
k_point[0], k_point[1], k_point[2]), "set xlabel 'T [K]'", "set ylabel 'E_{diff} [eV]'",
"plot 'bands_{}_{}_k_{}_{}_{}_difference_temperature.dat' u 1:2 w points pt 7 ps 1.2 lc rgb '#DC143C' notitle".format(lower_band_index,
lower_band_index+1, k_point[0], k_point[1], k_point[2])]

bd_change_plot_lines = ["set terminal postscript eps colour font 'Helvetica,20'", "set style data points",
"set output '| epstopdf --filter --outfile=bands_{}_{}_k_{}_{}_{}_difference_change_temperature.pdf'".format(lower_band_index,
lower_band_index+1, k_point[0], k_point[1], k_point[2]), "set xlabel 'T [K]'", "set ylabel '{/Symbol D}E_{diff} [meV]'",
"plot 'bands_{}_{}_k_{}_{}_{}_difference_temperature.dat' u 1:3 w points pt 7 ps 1.2 lc rgb '#DC143C' notitle".format(lower_band_index,
lower_band_index+1, k_point[0], k_point[1], k_point[2])]

# write to output
with open("bands_{}_{}_k_{}_{}_{}_difference_temperature.gnu".format(lower_band_index, lower_band_index+1, k_point[0], k_point[1],
                                                                     k_point[2]), "w") as plotfile:
    plotfile.write("\n".join(bd_plot_lines))

with open("bands_{}_{}_k_{}_{}_{}_difference_change_temperature.gnu".format(lower_band_index, lower_band_index+1, k_point[0], k_point[1],
                                                                            k_point[2]), "w") as plotfile:
    plotfile.write("\n".join(bd_change_plot_lines))
