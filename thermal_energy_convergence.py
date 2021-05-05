# file to write gnuplot scripts comparing the thermal energy calculated by Caesar at different supercell sizes
# written by Pascal Salzbrenner, pts28

# must be run in a directory containing Caesar calculations for different (at least one) supercell sizes, named supercell_xyz
# xyz indicate the supercell sizes
# each such supercell_xyz directory must contain an lte directory
# these contain the outputs of a Caesar harmonic lattice thermal energy calculation

import os
import sys

T_max = int(sys.argv[1])
num_atoms = int(sys.argv[2])

plotfile = open("thermal_energy_convergence.gnu", "w")

plotfile.write("set terminal postscript eps colour\n")
plotfile.write("set style data points\n")
plotfile.write("set key top left\n")
plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
plotfile.write("set output '| epstopdf --filter --outfile=thermal_energy_convergence.pdf'\n")

plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
plotfile.write("set linetype 2 lc rgb '#D95F02'\n")
plotfile.write("set linetype 3 lc rgb '#E6AB02'\n")
plotfile.write("set linetype 4 lc rgb '#66A61E'\n")
plotfile.write("set linetype 5 lc rgb '#191970'\n")
plotfile.write("set linetype 6 lc rgb '#7570B3'\n")
plotfile.write("set linetype 7 lc rgb '#E7298A'\n")
plotfile.write("set linetype 8 lc rgb '#1E90FF'\n")
plotfile.write("set linetype 9 lc rgb '#1B9E77'\n")
plotfile.write("set linetype 10 lc rgb '#B8860B'\n")
plotfile.write("set linetype cycle 10\n")

plotfile.write("set xlabel 'T [K]'\n")
plotfile.write("set ylabel 'F [meV]'\n")

plotfile.write("set xrange [0:{}]\n".format(T_max))

plot_string = "plot"

ls = os.listdir()

for dir in sorted(ls):
    if dir.startswith("supercell_"):
        # handle double-digit supercells
        supercell_size =  dir.split("_")
        if len(supercell_size) > 2:
            supercell_size_directions = supercell_size[1:]
        else:
            supercell_size_directions = [supercell_size[1][0], supercell_size[1][1], supercell_size[1][2]]

        plot_string += " '{}/lte/interpolated_free_energy.dat' u 1:($3*1000/{}) w points pt 7 ps 1.2 title '{}x{}x{} Supercell',".format(dir,
        num_atoms, supercell_size_directions[0], supercell_size_directions[1], supercell_size_directions[2])

plot_string = plot_string.rstrip(",")

plotfile.write("{}\n".format(plot_string))

plotfile.close()
