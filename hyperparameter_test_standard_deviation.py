# script to calculate the standard deviation of EDDP hyperparameter test MAE and RMSE values
# works with the hyperparameter_test.dat output from hyperparameter_test.sh

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

# open input and output files
infile = open("hyperparameter_test.dat", "r")
outfile = open("hyperparameter_test_standard_deviation.dat", "w")

# set up dictionary to contain the data
data_dict = {}

# read in the data
header = infile.readline()
outfile.write(header)

for line in infile:
    data = line.split()
    radial_cutoff = data[0]
    polynomials = data[1]
    rmse = float(data[3])
    mae = float(data[4])

    # add to the dictionary
    if radial_cutoff not in data_dict:
        data_dict[radial_cutoff] = {}

    if polynomials not in data_dict[radial_cutoff]:
        data_dict[radial_cutoff][polynomials] = [[], []]

    data_dict[radial_cutoff][polynomials][0].append(rmse)
    data_dict[radial_cutoff][polynomials][1].append(mae)

# calculate the standard deviation and write to the output file
for radial_cutoff in data_dict:
    for polynomials in data_dict[radial_cutoff]:
        rmse_std = np.std(data_dict[radial_cutoff][polynomials][0])
        mae_std = np.std(data_dict[radial_cutoff][polynomials][1])

        outfile.write("%s %s %f %f\n" % (radial_cutoff, polynomials, rmse_std, mae_std))

# close the files
infile.close()
outfile.close()
