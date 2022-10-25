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

allData = {}
for i,unRegDirectory in enumerate(IOlists['sfincsDir']):

    _, _, _, directory, _ = getFileInfo('/arbitrary/path', unRegDirectory, 'arbitrary')

    # Make target directory if it does not exist
    outDir = makeDir(saveDefaultTarget[i]) # Note that this script has file overwrite powers!
    
    # Retrieve the data
    dataFiles = findFiles('sfincsOutput.h5', directory) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in

    # Sort out the SFINCS subdirectories
    subdirs = [item.replace(directory+'/','') for item in dataFiles]

    splitSubdirs = [item.split('/') for item in subdirs]

    dataDepth = len(splitSubdirs[0])

    if not all([len(item) == dataDepth for item in splitSubdirs]):
        raise IOError('The structure of the SFINCS directory {} does not seem to be normal.'.format(directory))

    # Cycle through each file, read its data, and put that data in the proper place
    #radData = {}
    radDirName = None
    loadedData = {}
    for j,file in enumerate(dataFiles): # Scans through radial directories, and Er directories if present

        subdictNames = splitSubdirs[j]
        
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

        # Put the data in the appropriate place
        if dataDepth == 2: # Only radial directories are present
            allData[subdictNames[0]] = loadedData
        
        elif dataDepth == 3: # Radial and Er directories are present

            if radDirName == subdictNames[0]: # You are in the same radial directory as in the last iteration over dataFiles
                radData[subdictNames[1]] = loadedData
            
            else: # You are in a different radial directory from the last iteration over dataFiles
                if j != 0: # We shouldn't try to append data on the first loop iteration
                    allData[radDirName] = radData
                radDirName = subdictNames[0]
                radData = {}
                radData[subdictNames[1]] = loadedData
        
        else:
            raise IOError('The structure of the SFINCS directory {} does not seem to be normal.'.format(directory))

        loadedData = {} # This should be clean for each new file

    allData = {} # This should be clean for each new directory
