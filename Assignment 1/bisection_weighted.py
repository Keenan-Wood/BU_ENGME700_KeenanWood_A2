import numpy as np
from itertools import product

f = lambda x: (x-.101)**3
xBounds = np.array([-.1,1])

TOL = 10**-10
nIte_max = 10**5
xBuffer = .02

pt_coords = xBounds
pt_vals = np.zeros(2)
for i_pt in range(0,2):
    pt_vals[i_pt] = f(pt_coords[i_pt])

# Sort Coordinates and Values by value (ascending)
pt_coords = [val for _, val in sorted(zip(pt_vals, pt_coords))]
pt_vals.sort()

# Append midpoint Coordinate and Value
pt_coords = np.concatenate((pt_coords, [pt_coords[1]]))
pt_vals = np.concatenate((pt_vals, [pt_vals[1]]))

# Calculate coordinate of zero if function were linear
# -- Used to partition interval, instead of midpoint
# y(x) = slope * (x-x_0) + y(x_0)
# y(x) = (pt_vals[1]-pt_vals[0])/(pt_coords[1]-pt_coords[0]) * (x - pt_coords[0]) + pt_vals[0]
# calculate x so that y(x) = 0

nIte = 0
while nIte < nIte_max:
    nIte = nIte + 1
    pt_coords[2] = pt_coords[0] + min(1-xBuffer, max(xBuffer, -pt_vals[0]/ (pt_vals[1]-pt_vals[0]))) * (pt_coords[1]-pt_coords[0])
    pt_vals[2] = f(pt_coords[2])
    if abs(pt_vals[2]) <= TOL:
        break;
    elif pt_vals[2] < 0:
        pt_coords[0] = pt_coords[2]
        pt_vals[0] = pt_vals[2]
    else:
        pt_coords[1] = pt_coords[2]
        pt_vals[0] = pt_vals[2]

print("x=", pt_coords[2], ";", "nIte=", nIte)


TOL = 1