from __future__ import print_function
import os
import sys

try:
    from paraview.simple import OpenDataFile, RenameSource, Show
    pvloaded = True
except ImportError:
    pvloaded = False

for f in ["turbine", "walls"]:
    wd = "rvat-log-walls/" + f
    folders = sorted(os.listdir(wd))

    for folder in folders:
        files = os.listdir(wd + "/" + folder)
        files = sorted([int((file.replace("step_", "").split(".")[0])) for file in files])
        files = ["step_"+str(file)+".vtk" for file in files]
        if pvloaded:
            OpenDataFile([wd+"/"+folder+"/"+file for file in files])
            RenameSource(folder)
        else: print(files)
        
