# rescale the positions along the plotting k-path
# there is a mismatch between those from QE and from wannier90, this can rescale those from the former to match the latter
# it goes almost without saying that this program presupposes that the HSP order is equal in both cases
# requires as input bandstructures plotted using wannier90 and the QE bandplot.x utility

# naming conventions:
# QE and w90 calculations use the same fileroot
# the input for the pw.x bandstructure calculation must be named {fileroot}.bands
# gnuplot-plottable datafile from bandplot.x must be named {fileroot}_bands.dat.gnu

# written by Pascal Salzbrenner, pts28

import numpy as np

fileroot = input("What is the fileroot? ")

# open iput files
band_file = open("{}_bands.dat.gnu".format(fileroot), "r")
w90_hsp_file = open("{}_band.labelinfo.dat".format(fileroot), "r")
outfile = open("{}_QE_bands_rescaled.dat".format(fileroot), "w")

# read the number of points between HSPs in the QE calculation

qe_infile = open("{}.bands".format(fileroot), "r")
num_hsp_points = []

for line in qe_infile:

    if line.lstrip().lower().startswith("k_points"):
        # next line contains the number of HSPs
        num_hsp = int(qe_infile.readline())

        # the following num_hsp lines contain the HSPs - we are interested in the last element (assuming no comments) - the number of points along the lines
        for i in range(num_hsp):
            num_hsp_points.append(int(qe_infile.readline().split()[-1]))
    else:
        continue

qe_infile.close()

# read the wannier90 labelinfo file to determine the HSP locations in the FP band structure

w90_hsp_list = []

for line in w90_hsp_file:
    w90_hsp_list.append(float(line.split()[2]))

w90_hsp_file.close()

# we now read the {fileroot}_bands.dat.gnu
# we need two counters - one for which HSP we are at, and one for the number of points since the last HSP
hsp_counter = 0
points_since_hsp = 0

# we also a list of tuples containing (position, energy) data
read_data = []

for line in band_file:
    
    split_line = line.split()

    if split_line:

        data = np.array(split_line, dtype=float)

        if points_since_hsp == num_hsp_points[hsp_counter]:
            # we have read as many points since the last HSP as there are between it and the next HSP; ergo, the current point is an HSP
            # determine the distances between the two HSPs in the two different programs
            w90_hsp_dist = w90_hsp_list[hsp_counter+1] - w90_hsp_list[hsp_counter]
            qe_hsp_dist = data[0] - read_data[0][0]
            rescale_factor = w90_hsp_dist/qe_hsp_dist

            # write rescaled data to output
            for el in read_data:
                outfile.write("{} {}\n".format(el[0]*rescale_factor, el[1]))

            # the current HSP - we impose as its position the exact position from wannier90
            outfile.write("{} {}\n".format(w90_hsp_list[hsp_counter+1], data[1]))

            # set counters and data containers for the next set of points
            hsp_counter += 1
            points_since_hsp = 1
            read_data = [data]

        else:
            # regular point
            read_data.append(data)
            points_since_hsp += 1
    else:
        # the current line is empty - leave an empty line in the output as well
        outfile.write("\n")
        
        # this indicates that start of the block for the next eigenvalue - reset counters and containers to initial state
        hsp_counter = 0
        points_since_hsp = 0
        read_data = []

band_file.close()
outfile.close()
