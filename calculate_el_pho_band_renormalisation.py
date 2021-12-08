# postprocessing script to calculate the band gap renormalisation predicted by Caesar as a function of temperature
# written by Pascal Salzbrenner, pts28

# must be run in a directory containing a bs_valence_band and a bs_conduction_band directory, which contain the outputs of a Caesar bs
# calculation at the band edges for the valence and conduction bands respectively
# assumes both have the corrections calculated at the same temperatures

valence_band_renorm = open("bs_valence_band/band_gap_correction.dat", "r")
conduction_band_renorm = open("bs_conduction_band/band_gap_correction.dat", "r")

outfile = open("band_gap_renormalisation.dat", "w")

outfile.write("#Temperature [K], Change in band gap compared to static result [meV]\n")

for valence_line, conduction_line in zip(valence_band_renorm, conduction_band_renorm):
    outfile.write("{} {}\n".format(valence_line.split()[0], 1000*float(conduction_line.split()[1])-1000*float(valence_line.split()[1])))

valence_band_renorm.close()
conduction_band_renorm.close()
