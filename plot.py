# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, radialVarDict, adjustInputLengths, getFileInfo, makeDir, findFiles

# Set hard-coded reference variables
e = 1.602176634e-19 # C (proton charge)
#mBar = 1.67353e-27 # kg (hydrogen atom mass)
mBar = 1.672621911e-27 # kg (proton mass)
BBar = 1 # T
RBar = 1 # m
nBar = 1e20 # m^-3
TBar = 1.60217733e-16 # J = 1 keV
phiBar = 1000 # V = 1 kV
vBar = np.sqrt(2 * TBar / mBar) # m/s

# Get command line arguments and radial variables
args = getPlotArgs()
radialVars = radialVarDict()
radialVar = radialVars[args.radialVar[0]]
minBound = args.radialVarBounds[0]
maxBound = args.radialVarBounds[1]

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
    radDirName = None
    loadedData = {}
    for j,file in enumerate(dataFiles): # Scans through radial directories, and Er directories if present

        subdictNames = splitSubdirs[j]
        
        # Open the output file
        try:
            f = h5py.File(file, 'r')
        except:
            raise IOError('Unable to open SFINCS output (*.h5) file.')

        # Do a basic (not 100% conclusive) convergence check before reading an output file's data
        didNotConverge = []
        try:
            loadedData['finished'] = f['finished'][()]
        except KeyError:
            didNotConverge.append(file)
            continue

        # Read the desired data from the file
        defaults = ['Delta', 'alpha', 'nu_n']
        IVs = list(radialVars.values())
        DVs = ['Er', 'FSABjHat', 'FSABFlow', 'particleFlux_vm_psiHat', 'particleFlux_vm_psiN', 'particleFlux_vm_rHat', 'particleFlux_vm_rN', 'heatFlux_vm_psiHat',
                    'heatFlux_vm_psiN', 'heatFlux_vm_rHat', 'heatFlux_vm_rN', 'momentumFlux_vm_psiHat', 'momentumFlux_vm_psiN', 'momentumFlux_vm_rHat', 'momentumFlux_vm_rN']

        for varName in defaults + IVs + DVs:
            loadedData[varName] = f[varName][()]

        # Check that the default parameters are in order
        if loadedData['Delta'] != 0.0045694 or loadedData['alpha'] != 1 or loadedData['nu_n'] != 0.00833:
            raise IOError('It appears that the values of Delta, alpha, or nu_n were changed from their defaults. Please use the defaults to make unit conversions simpler.')

        # Put the data in the appropriate place
        if dataDepth == 2: # Only radial directories are present
            allData[subdictNames[0]] = loadedData
        
        elif dataDepth == 3: # Radial and Er directories are present

            if radDirName != subdictNames[0]: # You are in a different radial directory from the last iteration over dataFiles
                if j != 0: # We shouldn't try to append data on the first loop iteration
                    allData[radDirName] = radData
                radDirName = subdictNames[0]
                radData = {}
            
            radData[subdictNames[1]] = loadedData
        
        else:
            raise IOError('The structure of the SFINCS directory {} does not seem to be normal.'.format(directory))

        loadedData = {} # This should be clean for each new file

    # Now sort out what to plot
    stuffToPlot = {}
    for key,val in allData.items():
        
        if dataDepth == 2: # Only radial directories are present
            radialVal = val[radialVar]
        else: # Radial and Er directories are present
            radialVal = val[list(val.keys())[0]][radialVar] # Note that the same flux surface is used for each electric field sub-run
        
        minPass = minBound < 0 or minBound <= radialVal
        maxPass = maxBound < 0 or maxBound >= radialVal

        if minPass and maxPass:
            stuffToPlot[key] = val

    # Actually plot things
    # FIXME

     
    print(stuffToPlot)    
    allData = {} # This should be clean for each new directory

# Notify the user of convergence issues if necessary
if len(didNotConverge) > 0:
    print('It appears that the SFINCS run(s) which created the output file(s) in the following list did not converge properly:')
    print(didNotConverge)
