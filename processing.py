from __future__ import print_function
import os
import sys
import vtk
import numpy as np
import matplotlib.pyplot as plt
import timeseries as ts
from styleplot import styleplot

def load_vtk(t):
    f = "rvat-log/velocity/"+str(t)+".vtk"
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(f)
    reader.Update()
    data = reader.GetOutput()
    
    npoints = data.GetNumberOfPoints()
    point = data.GetPoint(0)
    d = data.GetPointData()
    
    array = d.GetArray("Velocity")
    
    u = np.zeros(npoints)
    v = np.zeros(npoints)
    w = np.zeros(npoints)
    x = np.zeros(npoints)
    y = np.zeros(npoints)
    z = np.zeros(npoints)
    
    for n in range(npoints):
        x[n], y[n], z[n] = data.GetPoint(n)
        u[n], v[n], w[n] = array.GetTuple(n)
        
    u = u[np.where(x==1.0)[0]]
    v = v[np.where(x==1.0)[0]]
    w = w[np.where(x==1.0)[0]]
    y = y[np.where(x==1.0)[0]]
    z = z[np.where(x==1.0)[0]]
        
    yarray, zarray, [uarray, varray, warray] = ts.build_plane_arrays(y, z,
                                                                    [u, v, w])
    return yarray, zarray, uarray, varray, warray
    
def calcwake():
    files = os.listdir("rvat-log/velocity")
    times = []
    for f in files:
        t = f.replace(".vtk", "")
        try:
            times.append(int(t))
        except ValueError:
            times.append(float(t))
    times.sort()
    t1 = 1.5
    i1 = times.index(t1)
    y, z, meanu, meanv, meanw = load_vtk(t1)
    i = 1
    for t in times[i1+1:]:
        i += 1
        y, z, u, v, w = load_vtk(t)
        meanu += u
        meanv += v
        meanw += w
    meanu = meanu/i
    meanv = meanv/i
    meanw = meanw/i
    return y, z, meanu, meanv, meanw
    
def plotwake():
    # Plot contours of mean streamwise velocity
    plt.figure(figsize=(10,5))
    y, z, u, v, w = calcwake()
    cs = plt.contourf(y/0.5, z, u, 20, cmap=plt.cm.coolwarm)
    plt.xlabel(r'$y/R$')
    plt.ylabel(r'$z/H$')
    styleplot()
    cb = plt.colorbar(cs, shrink=1, extend='both', 
                      orientation='horizontal', pad=0.3)
    cb.set_label(r'$U/U_{\infty}$')
    #turb_lines()
    ax = plt.axes()
    ax.set_aspect(2)
    plt.grid(True)
    plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
    plt.show()
    
def perf(plot=True):
    t, cp = np.loadtxt("rvat-log/performance.txt", unpack=True)
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
    print("Mean C_P =", np.mean(cp[i1:]))
    if plot:
        plt.figure()
        plt.plot(theta[ip1:], cp[ip1:])
        plt.xlabel(r"$\theta$ (deg)")
        plt.ylabel(r"$C_P$")
        styleplot()
        plt.show()

def main():
    plt.close("all")
    plotwake()
#    perf()

if __name__ == "__main__":
    main()
