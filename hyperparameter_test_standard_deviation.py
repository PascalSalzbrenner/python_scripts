# script to calculate the standard deviation of EDDP hyperparameter test MAE and RMSE values
# works with the hyperparameter_test.dat output from hyperparameter_test.sh

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

# open input and output files
infile = open("hyperparameter_test.dat", "r")

# set up dictionary to contain the data
data_dict = {}

# read in the data
# skip first line, which is just the header
infile.readline()

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

# close the input file
infile.close()

# calculate the standard deviation and write to the output file

outfile = open("hyperparameter_test_standard_deviation.dat", "w")

outfile.write("# -r [Ang]; -P; RMSE std [meV/atom]; MAE std [meV/atom]; RMSE [meV/atom]; MAE [meV/atom]\n")

for radial_cutoff in data_dict:
    for polynomials in data_dict[radial_cutoff]:
        rmse_std = np.std(data_dict[radial_cutoff][polynomials][0])
        mae_std = np.std(data_dict[radial_cutoff][polynomials][1])

        rmse_mean = np.mean(data_dict[radial_cutoff][polynomials][0])
        mae_mean = np.mean(data_dict[radial_cutoff][polynomials][1])

        outfile.write(f"{radial_cutoff} {polynomials} {rmse_std} {mae_std} {rmse_mean} {mae_mean}\n")

# close the output file
outfile.close()
