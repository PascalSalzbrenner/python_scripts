# it can happen in Caesar that the q-path lentghts for two different calculations of the same structure differ slightly
# eg for calculations using KS-DFT and OFDFT
# I think this is because the size of the unit cells is slightly different
# this script takes as input two directories in which there are completed Caesar calculations, and rescales the q-path of the second
# the longer k-path will be rescaled

# written by Pascal Salzbrenner, pts28

def clean_up_list_string(list_string):
    """Function to take a list of numbers and output a string containing them, separated by whitespace"""

    list_string_clean = ""

    for el in list_string:
        list_string_clean += "{} ".format(el)

    list_string_clean.rstrip()

    return list_string_clean



calculation_1 = input("What is the directory of the first calculation? ").rstrip("/")
calculation_2 = input("What is the directory of the second calculation? ").rstrip("/")
units = input("What units are the input dispersions in? ")

# determine where the HSPs are in the two different calculations
hsp_file_1 = open("{}/lte/phonon_dispersion.gnu".format(calculation_1), "r")
hsp_file_2 = open("{}/lte/phonon_dispersion.gnu".format(calculation_2), "r")

hsp_list_1 = []
hsp_list_2 = []

for line in hsp_file_1:
    if line.startswith("set xtics"):
        hsp_positions = line.split(",")
        break
    else:
        continue

for hsp_position in hsp_positions:
    hsp_position_raw = hsp_position.split()[-1]
    # only the last one will end in a bracket, but we can rstrip it for all of them with no consequences
    hsp_list_1.append(float(hsp_position_raw.rstrip(")")))

for line in hsp_file_2:
    if line.startswith("set xtics"):
        hsp_positions = line.split(",")
        break
    else:
        continue

for hsp_position in hsp_positions:
    hsp_position_raw = hsp_position.split()[-1]
    # only the last one will end in a bracket, but we can rstrip it for all of them with no consequences
    hsp_list_2.append(float(hsp_position_raw.rstrip(")")))

hsp_file_1.close()
hsp_file_2.close()

# determine which has the longer q-path:

if hsp_list_1[1] > hsp_list_2[1]:
    rescale_calc = calculation_1
    longer_list = hsp_list_1
    shorter_list = hsp_list_2
else:
    rescale_calc = calculation_2
    longer_list = hsp_list_2
    shorter_list = hsp_list_1

rescale_file = open("{}/lte/phonon_dispersion_curve_{}.dat".format(rescale_calc, units), "r")
rescale_outfile = open("{}/lte/phonon_dispersion_curve_{}_rescaled.dat".format(rescale_calc, units), "w")

# count how many HSPs we have come past
hsp_counter=1

# initialise the rescale factor
rescale_factor = shorter_list[1]/longer_list[1] # shorter_list[0]=longer_list[0] = 0

for line in rescale_file:
    position = float(line.split()[0])

    if position < longer_list[hsp_counter]:
        # non-HSP
        rescale_outfile.write("{} {}\n".format(shorter_list[hsp_counter-1]+(position-longer_list[hsp_counter-1])*rescale_factor,
        clean_up_list_string(line.split()[1:])))
    else:
        # HSP - enforce exact position
        rescale_outfile.write("{} {}\n".format(shorter_list[hsp_counter], clean_up_list_string(line.split()[1:])))

        # update counter and rescale factor
        if hsp_counter+1 < len(shorter_list):
            rescale_factor = (shorter_list[hsp_counter+1]-shorter_list[hsp_counter])/(longer_list[hsp_counter+1]-longer_list[hsp_counter])
            hsp_counter += 1
        else:
            # if this condition isn't fulfilled, we have reached the last HSP and hence the end of the file
            break

rescale_file.close()
rescale_outfile.close()
