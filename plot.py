# This script generates plots, plot data, and informational *.txt files given SFINCS output (*.h5) files. It can also perform basic convergence checks on the output files.
# Currently, this script cannot create 3D plots.
#FIXME WORK IN PROGRESS! fix all comments and such when finished. also note that this doesn't do 3d plotting. Also make sure that your input args explain everything properly

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, radialVarDict, adjustInputLengths, getFileInfo, makeDir, findFiles, writeFile, prettyRadialVar, prettyDataLabel, messagePrinter
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
    saveDefaultTarget = [join(item,'processed') for item in IOlists['sfincsDir']] # Give plots a subdirectory if no save locations are explicitely specified
else:
    saveDefaultTarget = IOlists['saveLoc'] 

# Create small functions that are useful only in this script
makeNeoclassicalNames = lambda x: [x+'_vm_'+IV for IV in IVs]
makeClassicalNames = lambda x: [x+'_'+IV for IV in IVs]
def writeInfoFile(listOfStrings, inputDir, outputDir, fileIDName):
    stringToWrite = ''.join(listOfStrings)
    fileToMake = join(outputDir, '{}-{}.txt'.format(inputDir, fileIDName))
    writeFile(fileToMake, stringToWrite, silent=True)

# Name some important variables
defaults = ['Delta', 'alpha', 'nu_n']

IVs = list(radialVars.values())
notRadialFluxes = ['Er', 'FSABjHat', 'FSABFlow']

particleFluxes = makeNeoclassicalNames('particleFlux')
heatFluxes = makeNeoclassicalNames('heatFlux')
momentumFluxes = makeNeoclassicalNames('momentumFlux')

classicalParticleFluxes = makeClassicalNames('classicalParticleFlux')
classicalParticleFluxesNoPhi1 = makeClassicalNames('classicalParticleFluxNoPhi1')
classicalHeatFluxes = makeClassicalNames('classicalHeatFlux')
classicalHeatFluxesNoPhi1 = makeClassicalNames('classicalHeatFluxNoPhi1')

nonCalcDVs = notRadialFluxes + particleFluxes + heatFluxes + momentumFluxes + classicalParticleFluxes + classicalParticleFluxesNoPhi1 + classicalHeatFluxes + classicalHeatFluxesNoPhi1
extras = ['Zs']

# Name some other variables to be calculated later (just radial current at the moment)
radialCurrents = []
for flux in particleFluxes:
   parts = flux.split('_')
   name = 'radialCurrent' + '_' + parts[1] + '_' + parts[2]
   radialCurrents.append(name)

DVs = nonCalcDVs + radialCurrents

# Loop through each directory
allData = {}
didNotConvergeAll = []
didNotConvergeDir = []
for i,unRegDirectory in enumerate(IOlists['sfincsDir']):
    
    # Regularize input directory name
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

        # Open the output file and do a basic (not 100% conclusive) convergence check before reading its data
        dirOfFileName = dirname(file)

        try:
            f = h5py.File(file, 'r')
            _ = f['finished'][()]
            shouldBePresent = f['FSABFlow'][()]
            if any(np.isnan(shouldBePresent)):
                raise ValueError
            convergenceState = 'PASS'
        
        except (IOError, KeyError, ValueError):
            didNotConvergeAll.append(file)
            didNotConvergeDir.append(file)
            convergenceState = 'FAIL'

        fileToWrite = join(dirOfFileName, 'convergence{}.txt'.format(convergenceState))
        convergenceString = 'This run {}ED basic convergence tests.'.format(convergenceState)
        writeFile(fileToWrite, convergenceString, silent=True)
        
        if args.checkConv or convergenceState == 'FAIL':
            continue

        # Read the desired data from the file
        for varName in defaults + IVs + nonCalcDVs + extras:
            loadedData[varName] = f[varName][()]

        # Check that the default parameters are in order
        if loadedData['Delta'] != 0.0045694 or loadedData['alpha'] != 1 or loadedData['nu_n'] != 0.00833:
            raise IOError('It appears that the values of Delta, alpha, or nu_n were changed from their defaults. Please use the defaults to make unit conversions simpler.')

        # Calculate other desired quantities (just radial current at the moment)
        for ind,radialCurrent in enumerate(radialCurrents):
           loadedData[radialCurrent] = np.dot(loadedData['Zs'], loadedData[particleFluxes[ind]])

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

    if not args.checkConv and len(didNotConvergeDir) != len(dataFiles):
    
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

        ErChoices = []
        IVvec = []
        for IV in IVs: # Select the radial variable you're plotting against
            
            DVvec = []
            for DV in DVs: # Select the data you want to plot
                
                baseName = nameOfDir + '-' + DV + '-vs-' + IV 
                plotName = baseName + '.pdf'
                dataName = baseName + '.dat'
                Zsname = baseName + '.Zs'

                fullPlotPath = join(outDir, plotName)
                fullDataPath = join(outDir, dataName)
                fullZsPath = join(outDir, Zsname)
        
                for radKey,radData in stuffToPlot.items():

                    if dataDepth == 2: # Only radial directories are present
                        dataToUse = radData

                    else: # Radial and Er directories are present
                        convergedErAbsVals = dict([(ErKey, np.abs(ErData['Er'])) for ErKey, ErData in radData.items()])
                        minErKey = min(convergedErAbsVals) # Returns key of Er subdirectory that has the smallest |Er| 
                        # If there are multiple minima, only the key of the first minimum will be returned.
                        # This should be fine - one would expect SFINCS runs with matching |Er| values to have converged to the same answer.
                        dataToUse = radData[minErKey]
                        ErChoices.append(join(radKey, minErKey) + '\n')
                   
                    IVvec.append(dataToUse[IV])
                    DVvec.append(fixOutputUnits(DV, dataToUse[DV]))

                IVvec = np.array(IVvec)
                DVvec = np.array(DVvec)
                
                DVshape = DVvec.shape
                if DVshape[-1] == 1: # Indicates that floats are being stores as single-element lists
                    DVvec = DVvec.reshape(DVshape[:-1]) # Gets rid of those extra lists so floats behave like floats

                combined = np.column_stack((IVvec,DVvec)) # The IV values will be the first column. The data comes in subsequent columns.
                combined = combined[combined[:, 0].argsort()] # This sorts the data so that radVar values are strictly ascending

                np.savetxt(fullDataPath, combined)

                plt.figure()
                plt.plot(combined[:,0], combined[:,1:]) # One horizontal axis data vector, (possibly) multiple vertical axis data vectors
                plt.xlabel(prettyRadialVar(IV))
                plt.ylabel(prettyDataLabel(DV))
                
                numLines = combined.shape[1] - 1
                
                if numLines > 1:
                    
                    Zs = dataToUse['Zs'] # Note that this assumes the Z for each species is the same throughout the plasma (i.e. the amount of stripping is constant)

                    np.savetxt(fullZsPath, Zs)
                    
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
        
        if len(didNotConvergeDir) > 0: # Note that if every output in an input directory did not converge, this file will not be written
            formattedList = [item + '\n' for item in didNotConvergeDir]
            writeInfoFile(formattedList, nameOfDir, outDir, 'didNotConverge')

        if len(ErChoices) > 0:
            uniqueChoices = list(set(ErChoices))
            uniqueChoices.sort()
            writeInfoFile(uniqueChoices, nameOfDir, outDir, 'ErChoices')
        
        allData = {} # This should be clean for each new directory
        didNotConvergeDir = [] # This should be clean for each new directory
        messagePrinter('Finished processing all available data in {}.'.format(directory))

# Notify the user of convergence issues if necessary
if len(didNotConvergeAll) > 0:
    messagePrinter('It appears that the SFINCS run(s) which created the output file(s) in the list below did not complete/converge properly.')
    print(didNotConvergeAll)
