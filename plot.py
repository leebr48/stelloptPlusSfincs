# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import h5py
import numpy as np

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, adjustInputLengths, getFileInfo, makeDir, findFiles

# Get command line arguments
args = getPlotArgs()

# Organize the directories that we will work in
inLists = {'sfincsDir':args.sfincsDir, 'saveLoc':args.saveLoc}
IOlists, _, _ = adjustInputLengths(inLists)

# If saveLoc was not specified, decide whether to use the profiles or equilibria locations
if all([item == None for item in IOlists['saveLoc']]):
    saveDefaultTarget = IOlists['sfincsDir']
else:
    saveDefaultTarget = IOlists['saveLoc']

for i,unRegDirectory in enumerate(IOlists['sfincsDir']):

    _, _, _, directory, _ = getFileInfo('/arbitrary/path', unRegDirectory, 'arbitrary')

    # Make target directory if it does not exist
    outDir = makeDir(actualSaveLoc) # Note that this script has file overwrite powers!

    # Retrieve the data
    dataFiles = findFiles('sfincsOutput.h5', directory) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in

    # Sort the data
    loadedData = {}
    ErData = {}
    radData = {} #FIXME new
    allData = {}
    for file in dataFiles: # Scans through radial directories, and Er directories if present
        
        # Open the output file
        try:
            f = h5py.File(file, 'r')
        except:
            raise IOError('Unable to open SFINCS output (*.h5) file.')

        # Read the desired data from the file
        loadedData['Er'] = f['Er'][()]
        loadedData['FSABjHat'] = f['FSABjHat'][()]
        loadedData['FSABFlow'] = f['FSABFlow'][()]
        loadedData['particleFlux_vm_psiN'] = f['particleFlux_vm_psiN'][()]
        loadedData['particleFlux_vm_rN'] = f['particleFlux_vm_rN'][()]
        loadedData['heatFlux_vm_psiN'] = f['heatFlux_vm_psiN'][()]
        loadedData['heatFlux_vm_rN'] = f['heatFlux_vm_rN'][()]
        loadedData['momentumFlux_vm_psiN'] = f['momentumFlux_vm_psiN'][()]
        loadedData['momentumFlux_vm_rN'] = f['momentumFlux_vm_rN'][()]

        # Put the data in the proper place
        fileNameToUse = file
        ErDirsIncluded = False

        while True:
            parentDirOfFile = dirname(fileNameToUse) # Walk up the file directory

            if parentDirOfFile == directory:
                # FIXME maybe now append the r bits?
                break # We want to work within the given directory, so we can stop going up the tree at this point
            
            nameOfParentDirOfFile = basename(parentDirOfFile)
            dirTypeIndicator = nameOfParentDirOfFile[0].lower()

            if dirTypeIndicator == 'p' or dirTypeIndicator == 'r': # You're in a radial directory
                
                if ErDirsIncluded:
                    # FIXME this appears to be busted... only one Er dictionary (the last one read) is ending up in each radial dictionary
                    # FIXME you need to add the ErData to another dictionary for the radius, and then append that dictionary to allData!
                    # FIXME can you just double-name dictionaries? like {{}}?
                    allData[nameOfParentDirOfFile] = ErData #FIXME we can't do this yet... we have to get the radial dic sorted first, THEN add it to allData
                    ErData = {}
                else:
                    allData[nameOfParentDirOfFile] = loadedData
            
            if dirTypeIndicator == 'e' or dirTypeIndicator == 'd': # You're in an Er directory
                ErDirsIncluded = True
                ErData[nameOfParentDirOfFile] = loadedData

            fileNameToUse = parentDirOfFile

        loadedData = {}
    print(allData)
