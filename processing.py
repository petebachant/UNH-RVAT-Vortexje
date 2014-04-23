from __future__ import print_function
import os
import sys
import vtk
import numpy as np
import matplotlib.pyplot as plt
import timeseries as ts
from styleplot import styleplot
import fdiff

R = 0.5
H = 1.0
U = 1.0

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
    
def calcwake(t1=1.5):
    files = os.listdir("rvat-log/velocity")
    times = []
    for f in files:
        t = f.replace(".vtk", "")
        try:
            times.append(int(t))
        except ValueError:
            times.append(float(t))
    times.sort()
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
    
def plotwake(plotlist=["meanu"], t1=1.5, save=False, savepath="", savetype=".pdf"):
    # Plot contours of mean streamwise velocity
    y, z, u, v, w = calcwake(t1)
    y_R = y/R
    z_H = z/H
    def turb_lines():
        plt.hlines(0.5, -1, 1, linestyles='solid', linewidth=2)
        plt.vlines(-1, 0, 0.5, linestyles='solid', linewidth=2)
        plt.vlines(1, 0, 0.5, linestyles='solid', linewidth=2)
    if "meanu" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
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
    if "meanv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y/0.5, z, v, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.3)
        cb.set_label(r'$V/U_{\infty}$')
        #turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.grid(True)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H, v, w, angles='xy')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.2, 0.1, r'$0.1$ m/s',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.tight_layout()
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'v-wquiver'+savetype)
    if "xvorticity" in plotlist or "all" in plotlist:
        dWdy = np.zeros(np.shape(u))
        dVdz = np.zeros(np.shape(u))
        for n in xrange(len(z)):
            dWdy[n,:] = fdiff.second_order_diff(w[n,:], y)
        for n in xrange(len(y)):
            dVdz[:,n] = fdiff.second_order_diff(v[:,n], z)
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, dWdy-dVdz, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.26)
        cb.set_label(r"$\Omega_x$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'xvorticity'+savetype)
    if "meancomboquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,6))
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, u, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cb.set_label(r'$U/U_{\infty}$')
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(y_R, z_H, v, w, angles='xy', width=0.0022)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.3, 0.1, r'$0.1 U_\infty$',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"meancomboquiv"+savetype)
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
    plotwake(plotlist=["meancomboquiv", "xvorticity"], t1=3)
#    perf()

if __name__ == "__main__":
    main()
