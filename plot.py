# This script generates plots, plot data, and informational *.txt files when fed SFINCS run directories. It can also perform basic convergence checks on the output files.
# Currently, this script cannot create 3D plots.
# To see the capabilities of this script, run it with the --help flag.

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, radialVarDict, adjustInputLengths, getFileInfo, makeDir, findFiles, writeFile, prettyRadialVar, prettyDataLabel, messagePrinter, now, saveTimeStampFile
from dataProc import checkConvergence, fixOutputUnits

# Get command line arguments and radial variables
args = getPlotArgs()
radialVars = radialVarDict()
radialVar = radialVars[args.radialVar[0]]
minBound = args.radialVarBounds[0]
maxBound = args.radialVarBounds[1]
IVs = list(radialVars.values())[:-1] # The last value is repeated (electric field definition)

# Organize the directories that we will work in
inLists = {'sfincsDir':args.sfincsDir, 'saveLoc':args.saveLoc}
IOlists, _, _ = adjustInputLengths(inLists)

# If saveLoc was not specified, decide whether to use the profiles or equilibria locations
if all([item == None for item in IOlists['saveLoc']]):
    saveDefaultTarget = [join(item,'processed') for item in IOlists['sfincsDir']] # Give plots a subdirectory if no save locations are explicitely specified
else:
    saveDefaultTarget = IOlists['saveLoc'] 

# Specify some small functions that are useful only in this script
makeNeoclassicalNames = lambda x: [x+distFunc+IV for IV in IVs]
makeOtherNames = lambda x: [x+'_'+IV for IV in IVs]
def writeInfoFile(listOfStrings, inputDir, outputDir, fileIDName):
    stringToWrite = ''.join(listOfStrings)
    fileToMake = join(outputDir, '{}-{}.txt'.format(inputDir, fileIDName))
    writeFile(fileToMake, stringToWrite, silent=True)
def findRadialInfo(radialDict, radialVar):
    try: # Will only work if there are no Er subdirectories
        radialVal = radialDict[radialVar]
        dataDepth = 2
    except KeyError:
        radialVal = radialDict[list(radialDict.keys())[0]][radialVar] # Note that the same flux surface is used in each Er subdirectory
        dataDepth = 3
    return radialVal, dataDepth

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
    dataFiles = findFiles('sfincsOutput.h5', directory, raiseError=True) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in

    # Cycle through each file, read its data, and put that data in the proper place
    radDirName = None
    loadedData = {}
    appendedSuccessfulData = 0
    for j,file in enumerate(dataFiles): # Scans through radial directories, and Er directories if present

        subdir = file.replace(directory+'/','')
        subdictNames = subdir.split('/')
        dataDepth = len(subdictNames)

        if dataDepth not in [2,3]:
            raise IOError('The structure of the SFINCS directory {} does not seem to be normal. This error was encountered while processing {}.'.format(directory, file))

        # Open the output file and do a basic (not 100% conclusive) convergence check before reading its data
        dirOfFileName = dirname(file)

        try:
            f = checkConvergence(file)
            convergenceState = 'PASS'
        
        except (IOError, KeyError):
            didNotConvergeAll.append(file)
            didNotConvergeDir.append(file)
            convergenceState = 'FAIL'

        convergenceStringList = ['File written ' + now() + '\n']
        convergenceStringList.append('This run {}ED basic convergence tests.\n'.format(convergenceState))
        writeInfoFile(convergenceStringList, basename(dirOfFileName), dirOfFileName, 'convergence')
        
        if args.checkConv:
            continue

        if convergenceState == 'PASS':
            # Check if we are in 'Phi1 mode' or not
            try:
                _ = f['Phi1Hat'][()]
                # Phi1 was included in the run
                distFunc = '_vd_'
            except KeyError:
                # Phi1 was not included in the run
                distFunc = '_vm_'

            # Name some important variables - doing this inside the loop is not particularly efficient, but it allows us to incorporate Phi1 effects automatically
            defaults = ['Delta', 'alpha', 'nu_n']

            notRadialFluxes = ['Er', 'FSABFlow', 'FSABjHat', 'FSABjHatOverRootFSAB2', 'FSABjHatOverB0']

            [neoclassicalParticleFluxes, neoclassicalHeatFluxes, neoclassicalMomentumFluxes] = [makeNeoclassicalNames(item) for item in ['particleFlux', 'heatFlux', 'momentumFlux']]

            [classicalParticleFluxes, classicalParticleFluxesNoPhi1, classicalHeatFluxes, classicalHeatFluxesNoPhi1] = [makeOtherNames(item) for item in ['classicalParticleFlux', 'classicalParticleFluxNoPhi1', 'classicalHeatFlux', 'classicalHeatFluxNoPhi1']]

            nonCalcDVs = notRadialFluxes + neoclassicalParticleFluxes + neoclassicalHeatFluxes + neoclassicalMomentumFluxes + classicalParticleFluxes + classicalParticleFluxesNoPhi1 + classicalHeatFluxes + classicalHeatFluxesNoPhi1

            extras = ['Zs', 'VPrimeHat']

            # Name some other variables to be calculated later
            [totalParticleFluxes, totalHeatFluxes] = [makeOtherNames(item) for item in ['totalParticleFlux', 'totalHeatFlux']]

            extensiveFluxes = ['extensiveParticleFlux', 'extensiveHeatFlux', 'extensiveMomentumFlux', 'extensiveClassicalParticleFlux', 'extensiveClassicalHeatFlux', 'extensiveTotalParticleFlux', 'extensiveTotalHeatFlux']

            radialCurrents = makeNeoclassicalNames('radialCurrent')
            extensiveRadialCurrent = ['extensiveRadialCurrent']

            DVs = nonCalcDVs + totalParticleFluxes + totalHeatFluxes + extensiveFluxes + radialCurrents + extensiveRadialCurrent
            
            # Read the desired data from the file
            for varName in defaults + IVs + extras: # These are all stored as scalars or 1D arrays that we want to keep as-is
                loadedData[varName] = f[varName][()] # The [()] converts the hdf5 object to a NumPy array

            for varName in nonCalcDVs:
                temp = f[varName][()]
                if varName == 'Er': # This is stored as a scalar
                    loadedData[varName] = temp
                else: # These are stored as 1D or 2D arrays
                    loadedData[varName] = temp[..., -1]

            # Check that the default parameters are in order
            if loadedData['Delta'] != 0.0045694 or loadedData['alpha'] != 1.0:
                raise IOError('It appears that the values of Delta or alpha were changed from their defaults. Please use the defaults to make unit conversions simpler.')

            # Calculate other desired quantities
            for radInd,(totalParticleFlux, totalHeatFlux) in enumerate(zip(totalParticleFluxes, totalHeatFluxes)):
                loadedData[totalParticleFlux] = loadedData[neoclassicalParticleFluxes[radInd]] + loadedData[classicalParticleFluxes[radInd]]
                loadedData[totalHeatFlux] = loadedData[neoclassicalHeatFluxes[radInd]] + loadedData[classicalHeatFluxes[radInd]]

            normalizedAreaFactor = loadedData['VPrimeHat'] # = dVHat/dpsiHat
            loadedData['extensiveParticleFlux'] = normalizedAreaFactor * loadedData['particleFlux'+distFunc+'psiHat']
            loadedData['extensiveHeatFlux'] = normalizedAreaFactor * loadedData['heatFlux'+distFunc+'psiHat']
            loadedData['extensiveMomentumFlux'] = normalizedAreaFactor * loadedData['momentumFlux'+distFunc+'psiHat']
            loadedData['extensiveClassicalParticleFlux'] = normalizedAreaFactor * loadedData['classicalParticleFlux_psiHat']
            loadedData['extensiveClassicalHeatFlux'] = normalizedAreaFactor * loadedData['classicalHeatFlux_psiHat']
            loadedData['extensiveTotalParticleFlux'] = loadedData['extensiveParticleFlux'] + loadedData['extensiveClassicalParticleFlux']
            loadedData['extensiveTotalHeatFlux'] = loadedData['extensiveHeatFlux'] + loadedData['extensiveClassicalHeatFlux']

            for radInd,radialCurrent in enumerate(radialCurrents):
                loadedData[radialCurrent] = np.dot(loadedData['Zs'], loadedData[neoclassicalParticleFluxes[radInd]])
            loadedData['extensiveRadialCurrent'] = normalizedAreaFactor * loadedData['radialCurrent'+distFunc+'psiHat']

            # Put the data in the appropriate place
            if radDirName != subdictNames[0]: # You are in a different radial directory from the last iteration over dataFiles
                if appendedSuccessfulData != 0: # We shouldn't try to append data before we've loaded it (not all variables are instantiated correctly yet)
                    allData[radDirName] = radData # Notice that in every loop iteration, we append data from the previous iteration

                radDirName = subdictNames[0]
                radData = {}
            
            if dataDepth == 2:
                radData = loadedData # With no Er directories, physical data is stored in the radial directories
            
            elif dataDepth == 3:
                radData[subdictNames[1]] = loadedData # With Er directories, physical data is stored inside of them, and they are inside the radial directories
            
            loadedData = {} # This should be clean for each new file
            appendedSuccessfulData += 1
            
        try:
            if j == len(dataFiles) - 1: # This is the last file, so we need to append the data before exiting the loop
                allData[radDirName] = radData

        except NameError: # This will only be entered if all the runs in a directory failed
            pass # The if statement directly below will handle this case
    
    if not args.checkConv and len(didNotConvergeDir) != len(dataFiles):
    
        # Now sort out what to plot
        stuffToPlot = {}
        for key,val in allData.items():
            
            radialVal, _ = findRadialInfo(val, radialVar)

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

                DVlower = DV.lower()
                if ('flux' in DVlower) and ('classical' not in DVlower) and ('total' not in DVlower): # If DV is a neoclassical flux
                    
                    nameEnd = DV[0].upper() + DV[1:]

                    if 'extensive' in DVlower:
                        DVnameForPlot = 'extensiveNeoclassical' + nameEnd.replace('Extensive','')
                    else:
                        DVnameForPlot = 'neoclassical' + nameEnd
                
                else:
                    DVnameForPlot = DV
                
                baseName = nameOfDir + '-' + DVnameForPlot + '-vs-' + IV 
                plotName = baseName + '.pdf'
                dataName = baseName + '.dat'
                Zsname = baseName + '.Zs'

                fullPlotPath = join(outDir, plotName)
                fullDataPath = join(outDir, dataName)
                fullZsPath = join(outDir, Zsname)
        
                for radKey,radData in stuffToPlot.items():
            
                    _, dataDepth = findRadialInfo(radData, IV)
                    
                    if dataDepth == 2: # Only radial directories are present
                        dataToUse = radData

                    else: # Radial and Er directories are present
                        convergedJrAbsVals = dict([(ErKey, np.abs(ErData['extensiveRadialCurrent'])[0]) for ErKey, ErData in radData.items()])
                        minJrKey = min(zip(convergedJrAbsVals.values(), convergedJrAbsVals.keys()))[1] # Returns key of Er subdirectory that has the smallest |Jr| 
                        dataToUse = radData[minJrKey]
                        ErChoices.append(join(radKey, minJrKey) + '\n')
                   
                    IVvec.append(dataToUse[IV])
                    DVvec.append(fixOutputUnits(DV, dataToUse[DV]))

                IVvec = np.array(IVvec)
                DVvec = np.array(DVvec)

                combined = np.column_stack((IVvec,DVvec)) # The IV values will be the first column. The data comes in subsequent columns.
                combined = combined[combined[:, 0].argsort()] # This sorts the data so that radVar values are strictly increasing

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

                if np.all(combined[:,0] > 0): # necessary because psiHat can be negative
                    plt.xlim(xmin=0)
                elif np.all(combined[:,0] < 0):
                    plt.xlim(xmax=0)
                
                plt.margins(0.01)
                
                plt.savefig(fullPlotPath, bbox_inches='tight', dpi=400)
                plt.close('all')

                IVvec = []
                DVvec = []
        
        if len(didNotConvergeDir) > 0: # Note that if every output in an input directory did not converge, this file will not be written
            formattedList = [item + '\n' for item in didNotConvergeDir]
            formattedList.insert(0, 'File written ' + now() + '\n')
            writeInfoFile(formattedList, nameOfDir, outDir, 'didNotConverge')

        if len(ErChoices) > 0:
            uniqueChoices = list(set(ErChoices))
            uniqueChoices.sort()
            uniqueChoices.insert(0, 'File written ' + now() + '\n')
            writeInfoFile(uniqueChoices, nameOfDir, outDir, 'ErChoices')
        
        allData = {} # This should be clean for each new directory
        didNotConvergeDir = [] # This should be clean for each new directory
        
        messagePrinter('Finished processing all available data in {}.'.format(directory))
        saveTimeStampFile(outDir, 'automatedPostprocessingLog', 'Data was last automatically postprocessed at this time: ')

# Notify the user of convergence issues if necessary
if len(didNotConvergeAll) > 0:
    messagePrinter('For your information: it appears that the SFINCS run(s) which created the output file(s) in the list below did not complete/converge.')
    messagePrinter(str(didNotConvergeAll))
