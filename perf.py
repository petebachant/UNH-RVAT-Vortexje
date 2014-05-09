#!/usr/bin/env python

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 0:
    case = sys.argv[1]
else:
    case = "free"

t, cp = np.loadtxt("rvat-log-"+case+"/performance.txt", unpack=True)

with open("rvat.cpp") as f:
    for line in f.readlines():
        line = line.split()
        try:
            if line[1] == "MILL_RADIUS":
                r = float(line[2])
            elif line[1] == "WIND_VELOCITY":
                uinfty = float(line[2])
            elif line[1] == "TIP_SPEED_RATIO":
                tsr = float(line[2])
        except IndexError:
            pass

omega = tsr*uinfty/r
theta = omega*t*180.0/np.pi
ip1 = 2
i1 = len(t)/2

plt.figure()
plt.plot(theta[ip1:], cp[ip1:])
plt.xlabel(r"$\theta$ (deg)")
plt.ylabel(r"$C_P$")
print("Mean C_P =", np.mean(cp[i1:]))
plt.show()
