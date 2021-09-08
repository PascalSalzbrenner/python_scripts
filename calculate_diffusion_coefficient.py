# script to calculate the diffusion coefficient based on a CASTEP AIMD run using ASE

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import sys
import ase.io, ase.md.analysis

# the user can specify how many timesteps are ignored, and into how many independent segments the trajectory is split for averaging
# these variables are directly passed to the calculate method of ase.md.analysis.DiffusionCoefficient
ignore_images = sys.argv[1]
num_segments = sys.argv[2]

# locate the .md file - this assumes that there is only one such file in the directory
ls = os.listdir()

# locate the relevant files
for file in ls:
    if file.endswith(".md"):
        md_file = file
    elif file.endswith(".param"):
        param_file = file

# read in every time step
trajectory = ase.io.read(md_file, index=":")

# set periodic boundary conditions
for i in range(len(trajectory)):
    trajectory[i].set_pbc(True)

# determine the atomic numbers and sort in ascending order
atomic_numbers = []
atomic_numbers.sort()

for atomic_number in trajectory[0].get_atomic_numbers():
    if atomic_number not in atomic_numbers:
        atomic_numbers.append(atomic_number)

# determine the timestep - assume it is given in fs - for other SI prefixes, an appropriate conversion factor can easily be defined
with open(param_file, "r") as infile:
    for line in infile:
        if "md_delta_t" in line.lower():
            castep_timestep = float(line.split()[2])

# define diffusion coefficient class, calculate and write out the diffusion coefficient
diffusion_coefficient = ase.md.analysis.DiffusionCoefficient(trajectory, timestep=castep_timestep*ase.units.fs)
diffusion_coefficient.calculate(ignore_n_images = ignore_images, number_of_segments = num_segments)

# this returns a list of lists
# each list is the diffusion coefficient for one element
# the first element of each element list is its diffusion coefficient, the second its standard deviation
# the elements are listed - I believe - according to atomic number
# the source code says "alphabetically", but I did a test with Yttrium and Lanthanum, where Yttrium is of course alphabetically after
# Lanthanum, but it has the lower atomic number, and Yttrium was returned first
diffusion_coefficients = diffusion_coefficient.get_diffusion_coefficients()

# convert to A^2/fs
for i in range(len(diffusion_coefficients)):
    for j in range(2):
        diffusion_coefficients[i][j] *= ase.units.fs

# write out diffusion coefficients
outfile = open("diffusion_coefficients",  "w")

outfile.write("Element number; Diffusion Coefficient [A^2/fs]; Standard Deviation [A^2/fs];  Diffusion Coefficient [cm^2/s]; Standard Deviation [cm^2/s]\n")

for i in range(len(diffusion_coefficients)):
    outfile.write("{} {} {} {} {}\n".format(atomic_numbers[i], diffusion_coefficients[i][0], diffusion_coefficients[i][1],
                                            diffusion_coefficients[i][0]*0.1, diffusion_coefficients[i][1]*0.1))

outfile.close()
