from __future__ import print_function
import os
import sys

try:
    from paraview.simple import OpenDataFile, RenameSource, Show
    pvloaded = True
except ImportError:
    pvloaded = False

casename = "free"

for f in ["turbine", "walls"]:
    wd = "rvat-log-" + casename + "/" + f
    folders = sorted(os.listdir(wd))
    for folder in folders:
        files = os.listdir(wd + "/" + folder)
        files = sorted([int((file.replace("step_", "").split(".")[0])) for file in files])
        files = ["step_"+str(file)+".vtk" for file in files]
        if pvloaded:
            OpenDataFile([wd+"/"+folder+"/"+file for file in files])
            RenameSource(folder)
        else: print(files)

# Load velocity
if os.path.isdir("rvat-log-" + casename + "/velocity"):
    wd = "rvat-log-" + casename + "/velocity"
    files = os.listdir(wd)
    filesn = []
    for f in files:
        try:
            filesn.append(int(f.replace(".vtk", "")))
        except ValueError:
            filesn.append(float(f.replace(".vtk", "")))
    files = [str(file)+".vtk" for file in sorted(filesn)]
    if pvloaded:
        OpenDataFile([wd+"/"+file for file in files])
        RenameSource("velocity")
    else: print(files)
