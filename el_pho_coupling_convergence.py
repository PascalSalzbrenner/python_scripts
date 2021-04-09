# file to write gnuplot scripts comparing the energies of the valence and conduction bands at different supercell sizes
# written by Pascal Salzbrenner, pts28

# must be run in a directory containing Caesar calculations for different (at least one) supercell sizes, named supercell_xyz
# xyz indicate the supercell sizes
# each such supercell_xyz directory must contain a bs_valence_band and a bs_conduction_band directory
# these contain the outputs of a Caesar bs calculation at the band edges for the valence and conduction bands respectively
# if band_gap_renormalisation.dat is also present, the band gap at different supercell sizes is calculated as well

import os

vb_plotfile = open("valence_band_renormalisation.gnu", "w")
cb_plotfile = open("conduction_band_renormalisation.gnu", "w")
bg_plotfile = open("band_gap_renormalisation.gnu", "w")

vb_plotfile.write("set terminal postscript eps colour\n")
vb_plotfile.write("set style data points\n")
vb_plotfile.write("set key top left\n")
vb_plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
vb_plotfile.write("set output '| epstopdf --filter --outfile=valence_band_renormalisation_convergence.pdf'\n")

cb_plotfile.write("set terminal postscript eps colour\n")
cb_plotfile.write("set style data points\n")
cb_plotfile.write("set key top left\n")
cb_plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
cb_plotfile.write("set output '| epstopdf --filter --outfile=conduction_band_renormalisation_convergence.pdf'\n")

bg_plotfile.write("set terminal postscript eps colour\n")
bg_plotfile.write("set style data points\n")
bg_plotfile.write("set key top left\n")
bg_plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
bg_plotfile.write("set output '| epstopdf --filter --outfile=band_gap_renormalisation_convergence.pdf'\n")

vb_plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
vb_plotfile.write("set linetype 2 lc rgb '#FF8811'\n")
vb_plotfile.write("set linetype 3 lc rgb '#FFCC00'\n")
vb_plotfile.write("set linetype 4 lc rgb '#66A61E'\n")
vb_plotfile.write("set linetype 5 lc rgb '#191970'\n")
vb_plotfile.write("set linetype 6 lc rgb '#7570B3'\n")
vb_plotfile.write("set linetype 7 lc rgb '#E7298A'\n")
vb_plotfile.write("set linetype 8 lc rgb '#1E90FF'\n")
vb_plotfile.write("set linetype 9 lc rgb '#1B9E77'\n")
vb_plotfile.write("set linetype 10 lc rgb '#B8860B'\n")
vb_plotfile.write("set linetype cycle 10\n")

cb_plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
cb_plotfile.write("set linetype 2 lc rgb '#FF8811'\n")
cb_plotfile.write("set linetype 3 lc rgb '#FFCC00'\n")
cb_plotfile.write("set linetype 4 lc rgb '#66A61E'\n")
cb_plotfile.write("set linetype 5 lc rgb '#191970'\n")
cb_plotfile.write("set linetype 6 lc rgb '#7570B3'\n")
cb_plotfile.write("set linetype 7 lc rgb '#E7298A'\n")
cb_plotfile.write("set linetype 8 lc rgb '#1E90FF'\n")
cb_plotfile.write("set linetype 9 lc rgb '#1B9E77'\n")
cb_plotfile.write("set linetype 10 lc rgb '#B8860B'\n")
cb_plotfile.write("set linetype cycle 10\n")

bg_plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
bg_plotfile.write("set linetype 2 lc rgb '#FF8811'\n")
bg_plotfile.write("set linetype 3 lc rgb '#FFCC00'\n")
bg_plotfile.write("set linetype 4 lc rgb '#66A61E'\n")
bg_plotfile.write("set linetype 5 lc rgb '#191970'\n")
bg_plotfile.write("set linetype 6 lc rgb '#7570B3'\n")
bg_plotfile.write("set linetype 7 lc rgb '#E7298A'\n")
bg_plotfile.write("set linetype 8 lc rgb '#1E90FF'\n")
bg_plotfile.write("set linetype 9 lc rgb '#1B9E77'\n")
bg_plotfile.write("set linetype 10 lc rgb '#B8860B'\n")
bg_plotfile.write("set linetype cycle 10\n")

vb_plotfile.write("set xlabel 'T [K]'\n")
vb_plotfile.write("set ylabel '{/Symbol D}E_{Band} [meV]'\n")

cb_plotfile.write("set xlabel 'T [K]'\n")
cb_plotfile.write("set ylabel '{/Symbol D}E_{Band} [meV]'\n")

bg_plotfile.write("set xlabel 'T [K]'\n")
bg_plotfile.write("set ylabel '{/Symbol D}E_{G} [meV]'\n")

vb_plot_string = "plot"
cb_plot_string = "plot"
bg_plot_string = "plot"

ls = os.listdir()

for dir in sorted(ls):
    if dir.startswith("supercell_"):
        # handle double-digit supercells
        supercell_size =  dir.split("_")
        if len(supercell_size) > 2:
            supercell_size_directions = supercell_size[1:]
        else:
            supercell_size_directions = [supercell_size[1][0], supercell_size[1][1], supercell_size[1][2]]

        if "bs_valence_band" in os.listdir(dir):
            vb_plot_string += " '{}/bs_valence_band/band_gap_correction.dat' u 1:($2*1000) w points pt 7 ps 1.2 title '{}x{}x{} Supercell',".format(dir,
            supercell_size_directions[0], supercell_size_directions[1], supercell_size_directions[2])
        if "bs_conduction_band" in os.listdir(dir):
            cb_plot_string += " '{}/bs_conduction_band/band_gap_correction.dat' u 1:($2*1000) w points pt 7 ps 1.2 title '{}x{}x{} Supercell',".format(dir,
            supercell_size_directions[0], supercell_size_directions[1], supercell_size_directions[2])
        if "band_gap_renormalisation.dat" in os.listdir(dir):
            bg_plot_string += " '{}/band_gap_renormalisation.dat' u 1:2 w points pt 7 ps 1.2 title '{}x{}x{} Supercell',".format(dir,
            supercell_size_directions[0], supercell_size_directions[1], supercell_size_directions[2])

vb_plot_string = vb_plot_string.rstrip(",")
cb_plot_string = cb_plot_string.rstrip(",")
bg_plot_string = bg_plot_string.rstrip(",")

vb_plotfile.write("{}\n".format(vb_plot_string))
cb_plotfile.write("{}\n".format(cb_plot_string))
bg_plotfile.write("{}\n".format(bg_plot_string))

vb_plotfile.close()
cb_plotfile.close()
bg_plotfile.close()
