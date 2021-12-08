# code to determine the convex hull enveloping a list of q-points
# in the limit of high grid density, this ought to be equivalent to the actual vertices of the IBZ

# input: paths to directories containing the desired Caesar-generated ibz.dat files

import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D

# set up q-point list
points = []

# read in q-points int he IBZ
ibz_file = open("ibz.dat", "r")

for line in ibz_file:
    split_line=line.split()
    points.append(np.array(split_line[:-1], dtype=np.float64))

# create convex hull
hull = ConvexHull(points)

# extract vectors describing vertices
vertices=[]

for index in hull.vertices:
    vertices.append(points[index])

# print vertices
print(len(hull.vertices))
print(vertices)

with open("ibz_vertices.txt", "w") as vertex_file:
    vertex_file.write("Number of vertices: {}\n".format(len(hull.vertices)))

    for vertex in vertices:
        write_string = "{}".format(vertex[0])

        for el in vertex[1:]:
            write_string += " {}".format(el)

        write_string += "\n"

        vertex_file.write(write_string)

vertices = np.array(vertices)
points = np.array(points)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot(points.T[0], points.T[1], points.T[2], "bo", zorder=1)

for s in hull.simplices:
    s = np.append(s, s[0])
    ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-", zorder=2
    )

ax.plot(vertices.T[0], vertices.T[1], vertices.T[2], color="red", marker="o", linestyle="")

plt.savefig("ibz.pdf")
plt.show()
