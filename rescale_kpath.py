# rescale the positions along the plotting k-path
# there is a mismatch between those from my own code and from VASP/wannier90, this can rescale those from the latter to match the former
# it goes almost without saying that this program presupposes that the HSP order is equal in both cases
# requires as input a bandstructure plotted using either sumo or wannier90
# written by Pascal Salzbrenner, pts28

import numpy as np

code = input("Was the first-principles band structure generated using sumo or wannier90? ").lower()

# open iput files
if code.startswith("s"):
    # generated using sumo
    band_file = open("band.dat", "r")
    fp_hsp_file = open("sumo-bandplot.log", "r")
    outfile = open("band_rescaled.dat", "w")
elif code.startswith("w"):
    # generated using wannier90
    band_file = open("wannier90_band.dat", "r")
    fp_hsp_file = open("wannier90_band.gnu", "r")
    outfile = open("wannier90_band_rescaled.dat", "w")
else:
    raise ValueError("Only sumo and wannier90 generated band structures are supported. Please supply one or the other and indicate which.")

# read the CamTB HSP positions

hsp_file = open("hsp.dat", "r")
camtb_hsp_list = []

for line in hsp_file:
    if not line.startswith("#"):
        camtb_hsp_list.append(np.double(line.split()[1]))

hsp_file.close()

# read the fp_hsp_file file to determine the HSP locations in the FP band structure
# how it must be read depends on which code was used to plot the FP band structure

fp_hsp_list = []

if code.startswith("s"):

    for i in range(2):
        # read past the first two lines
        fp_hsp_file.readline()

    # the rest of the lines in sumo-bandplot.log contain an HSP each

    for line in fp_hsp_file:
        if not line.split():
            # the line is empty (final line)
            continue

        fp_hsp_list.append(float(line.split(":")[0]))

elif code.startswith("w"):
    for line in fp_hsp_file:
        if line.startswith("set xtics"):
            hsp_positions = line.split(",")
            break
        else:
            continue

    for hsp_position in hsp_positions:
        hsp_position_raw = hsp_position.split()[-1]
        if hsp_position_raw.endswith(")"):
            # the last one will end with the bracket
            fp_hsp_list.append(np.double(hsp_position_raw.rstrip(")")))
        else:
            fp_hsp_list.append(np.double(hsp_position_raw))

fp_hsp_file.close()

for line in band_file:

    if line.startswith("#"):
        # first line in sumo-generated band.dat
        continue

    data = line.split()

    if data:
        isclose_array = np.isclose(np.double(data[0]), np.array(fp_hsp_list), atol=1e-5)
        if isclose_array.any():
            # the current point is one of the HSPs - isclose_array will have only a single True element, which gives the HSP index
            hsp_index = np.where(isclose_array==True)[0][0]
            prev_camtb_hsp = camtb_hsp_list[hsp_index]
            prev_w90_hsp = fp_hsp_list[hsp_index]
            outfile.write("{} {}\n".format(prev_camtb_hsp, data[1])) # impose that the HSP has the position
            # of the corresponding CamTB HSP
            if hsp_index != len(fp_hsp_list)-1:
                # not the last HSP
                dist_w90 = fp_hsp_list[hsp_index+1] - prev_w90_hsp
                dist_camtb = camtb_hsp_list[hsp_index+1] - prev_camtb_hsp
                rescale_factor = dist_camtb/dist_w90
        else:
            # the current point is one between HSPs
            outfile.write("{} {}\n".format(prev_camtb_hsp+rescale_factor*(np.double(data[0])-prev_w90_hsp), data[1]))
    else:
        outfile.write("\n")

band_file.close()
outfile.close()
