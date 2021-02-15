# python script to transform, given two bases X and Y and a number of vectors, the representation of those vectors from the first to the second basis

# input file: basis_transform.in
# format:
# X_11 X_12 X_13
# X_21 X_22 X_23
# X_31 X_32 X_33
#
# Y_11 Y_12 Y_13
# Y_21 Y_22 Y_23
# Y_31 Y_32 Y_33
#
# a_x1 a_x2 a_x3
# b_x1 b_x2 b_x3
# ...
# z_x1 z_x2 z_x3

# each row in X and Y is a basis vector
# only the first three (whitespace-separated) elements of every line will be read - the rest can be anything
# there must be no comment full comment lines anywhere

# output: basis_transform.out, containing the vectors
# a_y1 a_y2 a_y3
# b_y1 b_y2 b_y3
# ...
# z_y1 z_y2 Z_y3

# written by Pascal Salzbrenner pts28@cam.ac.uk

import numpy as np

# input and output files
infile = open("basis_transform.in", "r")
outfile = open("basis_transform.out", "w")

basis_X = np.array([infile.readline().split()[:3], infile.readline().split()[:3], infile.readline().split()[:3]], dtype=float)

# empty line
infile.readline()

basis_Y = np.array([infile.readline().split()[:3], infile.readline().split()[:3], infile.readline().split()[:3]], dtype=float)

basis_transform = np.matmul(np.linalg.inv(basis_Y.T), basis_X.T)

# empty line
infile.readline()

# list for the output vectors
vectors_Y = []

for line in infile:
    vector_X = np.array(line.split()[:3], dtype=float)
    vector_Y = np.dot(basis_transform, vector_X)
    outfile.write("{} {} {}\n".format(vector_Y[0], vector_Y[1], vector_Y[2]))

infile.close()
outfile.close()
