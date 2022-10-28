# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished. also note that this doesn't do 3d plotting 

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, radialVarDict, adjustInputLengths, getFileInfo, makeDir, findFiles, prettyRadialVar, prettyDataLabel
from dataProc import fixOutputUnits

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
didNotConverge = []
for i,unRegDirectory in enumerate(IOlists['sfincsDir']):

    _, _, _, directory, _ = getFileInfo('/arbitrary/path', unRegDirectory, 'arbitrary')

    # Make target directory if it does not exist
    outDir = makeDir(saveDefaultTarget[i]) # Note that this script has file overwrite powers!
    
    # Retrieve the data
    dataFiles = findFiles('sfincsOutput.h5', directory) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in

    if len(dataFiles) == 0:
        raise IOError('No SFINCS output (*.h5) file could be found in the input directory {}.'.format(directory))

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
        try:
            _ = f['finished'][()]
            shouldBePresent = f['FSABFlow'][()]
        except KeyError:
            didNotConverge.append(file)
            continue

        if any(np.isnan(shouldBePresent)):
            didNotConverge.append(file)
            continue

        if args.checkConv:
            continue

        # Read the desired data from the file
        defaults = ['Delta', 'alpha', 'nu_n']
        IVs = list(radialVars.values())
        DVs = ['Er', 'FSABjHat', 'FSABFlow', 'particleFlux_vm_psiHat', 'particleFlux_vm_psiN', 'particleFlux_vm_rHat', 'particleFlux_vm_rN', 'heatFlux_vm_psiHat',
                    'heatFlux_vm_psiN', 'heatFlux_vm_rHat', 'heatFlux_vm_rN', 'momentumFlux_vm_psiHat', 'momentumFlux_vm_psiN', 'momentumFlux_vm_rHat', 'momentumFlux_vm_rN']
        extras = ['Zs']

        for varName in defaults + IVs + DVs + extras:
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

    if not args.checkConv:
    
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
        nameOfDir = basename(directory)

        if dataDepth == 2: # Only radial directories are present #FIXME might need to put this logic elsewhere, hopefully in a more compressed/general way than repeating loops

            IVvec = []
            for IV in IVs: # Select the radial variable you're plotting against
                
                DVvec = []
                for DV in DVs: # Select the data you want to plot
                        
                    plotName = nameOfDir + '-' + DV + '-vs-' + IV + '.pdf'

                    fullPlotPath = join(outDir, plotName)
            
                    for key,data in stuffToPlot.items():
                       
                        IVvec.append(data[IV])
                        DVvec.append(fixOutputUnits(DV, data[DV]))

                    IVvec = np.array(IVvec)
                    DVvec = np.array(DVvec)
                    
                    DVshape = DVvec.shape
                    if DVshape[-1] == 1: # Indicates that floats are being stores as single-element lists
                        DVvec = DVvec.reshape(DVshape[:-1]) # Gets rid of those extra lists so floats behave like floats

                    combined = np.column_stack((IVvec,DVvec)) # The IV values will be the first column. The data comes in subsequent columns.
                    combined = combined[combined[:, 0].argsort()] # This sorts the data so that radVar values are strictly ascending

                    plt.figure()
                    plt.plot(combined[:,0], combined[:,1:]) # One horizontal axis data vector, (possibly) multiple vertical axis data vectors
                    plt.xlabel(prettyRadialVar(IV))
                    plt.ylabel(prettyDataLabel(DV))
                    
                    numLines = combined.shape[1] - 1
                    
                    if numLines > 1:
                        
                        Zs = stuffToPlot[list(stuffToPlot.keys())[0]]['Zs'] # Note that this assumes the Z for each species is the same throughout the plasma (i.e. the amount of stripping is constant)
                        
                        leg = []
                        for specNum in range(numLines):
                            leg.append(r'$Z={}$'.format(int(Zs[specNum])))

                        plt.legend(leg, loc='best')

                    plt.xlim(xmin=0)
                    plt.margins(0.01)
                    
                    plt.savefig(fullPlotPath, bbox_inches='tight', dpi=400)
                    plt.close('all')

                    IVvec = []
                    DVvec = []
                        
        allData = {} # This should be clean for each new directory

# Notify the user of convergence issues if necessary
if len(didNotConverge) > 0:
    print('It appears that the SFINCS run(s) which created the output file(s) in the following list did not complete/converge properly:')
    print(didNotConverge)
