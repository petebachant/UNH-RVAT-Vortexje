from __future__ import print_function
import os
import sys

try:
    from paraview.simple import OpenDataFile, RenameSource, Show
    pvloaded = True
except ImportError:
    pvloaded = False

try:
    casename = sys.argv[1]
except AttributeError:
    casename = "rvat"

wd = casename + "-log/" + casename
folders = sorted(os.listdir(wd))

for folder in folders:
    files = os.listdir(wd + "/" + folder)
    files = sorted([int((file.replace("step_", "").split(".")[0])) for file in files])
    files = ["step_"+str(file)+".vtk" for file in files]
    if pvloaded:
        OpenDataFile([wd+"/"+folder+"/"+file for file in files])
        RenameSource(folder)
    else: print(files)
    
