# file process the output of fourier_interpolation calculations with Caesar (bm418) where some high symmetry points in the phonon dispersion
# are split, and generate the regular output files generated from such a calculation when the HSPs aren't split
# written by Pascal Salzbrenner - pts28

# input files: path_*.dat, high_symmetry_points_*.dat, phonon_dispersion_curve_*.dat

import sys

def list_to_string(input_list):
    """Function that takes a list and generates a string that consists of the list's elements with a whitespace between them"""

    return_string = ""

    for el in input_list:
        return_string += "{} ".format(el)

    return return_string.rstrip()

# read the index of last file - the first file is always assumed to be 1
last_file = int(sys.argv[1])

# output files

outpath = open("path.dat", "w")
outhsp = open("high_symmetry_points.dat", "w")
outcurve = open("phonon_dispersion_curve.dat", "w")

# the positions and numbers in the files have to be shifted when they are concatenated, these variables determine the shift

add_hsp = 0
add_hsp_pos = 0
add_curve_pos = 0

# iterate over all files
for i in range(1, last_file+1):

    # at the start of each iteration, these counters are set to 0 so the program knows the line it is reading is the first
    path_count = 0
    hsp_count = 0
    curve_count = 0

    with open("path_{}.dat".format(i), "r") as pathfile:

        for line in pathfile:
            if path_count == 0:
                if i == 1:
                    # the first file
                    outpath.write(line)
                else:
                    # this is done for any file but the first
                    outpath.write("{} | {} # {}|{}\n".format(prev_path_line.split("#")[0].strip(), line.split("#")[0].strip(),
                                                           prev_path_line.split("#")[1].strip(), line.split("#")[1].strip()))
                prev_path_line = pathfile.readline()
            else:
                outpath.write(prev_path_line)
                prev_path_line = line
            path_count +=1

    with open("high_symmetry_points_{}.dat".format(i), "r") as hspfile:

        for line in hspfile:
            if hsp_count == 0:
                if i == 1:
                    outhsp.write(line)
                else:
                    outhsp.write("{} {}\n".format(int(prev_hsp_line.split()[0])+add_hsp, float(prev_hsp_line.split()[1])+add_hsp_pos))
                    # the previous line and this one are the same position in the graph so this line isn't written
                    # update the offsets
                    add_hsp += (int(prev_hsp_line.split()[0])-1) # because the first line of every file isn't written, we add one less to every hsp
                    # otherwise there would be numbers which are skipped
                    add_hsp_pos += float(prev_hsp_line.split()[1])
                prev_hsp_line = hspfile.readline()
            else:
                outhsp.write("{} {}\n".format(int(prev_hsp_line.split()[0])+add_hsp, float(prev_hsp_line.split()[1])+add_hsp_pos))
                prev_hsp_line = line
            hsp_count += 1

    with open("phonon_dispersion_curve_{}.dat".format(i), "r") as curvefile:

        for line in curvefile:
            if curve_count == 0:
                if i == 1:
                    outcurve.write(line)
                    prev_curve_line = curvefile.readline()
                else:
                    outcurve.write("{} {}\n".format(float(prev_curve_line.split()[0])+add_curve_pos, list_to_string(prev_curve_line.split()[1:])))
                    # update the offset
                    add_curve_pos += float(prev_curve_line.split()[0])
                    prev_curve_line = line
            else:
                outcurve.write("{} {}\n".format(float(prev_curve_line.split()[0])+add_curve_pos, list_to_string(prev_curve_line.split()[1:])))
                prev_curve_line = line
            curve_count += 1

outpath.write(prev_path_line)
outhsp.write("{} {}\n".format(int(prev_hsp_line.split()[0])+add_hsp, float(prev_hsp_line.split()[1])+add_hsp_pos))
outcurve.write("{} {}\n".format(float(prev_curve_line.split()[0])+add_curve_pos, list_to_string(prev_curve_line.split()[1:])))

outpath.close()
outhsp.close()
outcurve.close()
