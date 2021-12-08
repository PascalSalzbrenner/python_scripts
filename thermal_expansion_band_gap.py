# Python script which writes gnuplot files for plotting the (change in) band gap as a function of temperature

# input
# band_gap_temperature.dat file, with the format T [K], band gap [eV], change in band gap compared to static [meV]

# output
# band_gap_temperature.gnu, to plot the absolute value of the band gap as a function of temperature
# band_gap_change_temperature.gnu, to plot the change in the band gap compared to the static result as a function of temperature

# create lists of lines for band gap gnuplot files
bg_plot_lines = ["set terminal postscript eps colour font 'Helvetica,20'", "set style data points",
"set output '| epstopdf --filter --outfile=band_gap_temperature.pdf'", "set xlabel 'T [K]'", "set ylabel 'E_{G} [eV]'",
"plot 'band_gap_temperature.dat' u 1:2 w points pt 7 ps 1.2 lc rgb '#DC143C' notitle"]

bg_change_plot_lines = ["set terminal postscript eps colour font 'Helvetica,20'", "set style data points",
"set output '| epstopdf --filter --outfile=band_gap_change_temperature.pdf'", "set xlabel 'T [K]'", "set ylabel '{/Symbol D}E_{G} [meV]'",
"plot 'band_gap_temperature.dat' u 1:3 w points pt 7 ps 1.2 lc rgb '#DC143C' notitle"]

# write to output
with open("band_gap_temperature.gnu", "w") as plotfile:
    plotfile.write("\n".join(bg_plot_lines))

with open("band_gap_change_temperature.gnu", "w") as plotfile:
    plotfile.write("\n".join(bg_change_plot_lines))
