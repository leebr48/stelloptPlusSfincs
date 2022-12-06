# FIXME This script is used to copy input.namelist files and convert them to running with Phi1.
# FIXME Maybe mention that this script will pick the best Er run (if applicable) and only copy that one?
# FIXME if you don't copy over equilibrium and such, note that it/they cannot move!
# FIXME should you check that the loaded runs don't already have Phi1?

from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
import numpy as np
thisFile = abspath(getfile(currentframe()))
thisDir = dirname(thisFile)
sys.path.append(join(thisDir, 'src/'))
from IO import getPhi1SetupArgs, getFileInfo, adjustInputLengths, makeDir, findFiles, writeFile
from dataProc import checkConvergence
_, thisFileName, _, _, _ = getFileInfo(thisFile, 'arbitrary/path', 'arbitrary')

# Get command line arguments
args = getPhi1SetupArgs()

# Organize the directories that we will work in and get their paths if necessary
inDirs = [getFileInfo('/arbitrary/path', unRegDirectory, 'arbitrary')[3] for unRegDirectory in args.sfincsDir]
regularizedSaveLocs = [getFileInfo('/arbitrary/path', unRegDirectory, 'arbitrary')[3] if unRegDirectory is not None else None for unRegDirectory in args.saveLoc]
inLists = {'sfincsDir':inDirs, 'saveLoc':regularizedSaveLocs}
IOlists, _, _ = adjustInputLengths(inLists)

# If saveLoc was not specified, use default location
if all([item == None for item in IOlists['saveLoc']]):
    outDirs = [item + '_Phi1' for item in IOlists['sfincsDir']]
else:
    outDirs = IOlists['saveLoc'] 

# Set variables that need to be changed in input.namelist
newRunParams = {'ambipolarSolve': {'val': '.false.', 'paramList': 'general', 'used':False},
                'includePhi1': {'val': '.true.', 'paramList': 'physicsParameters', 'used':False}
                }

# I/O flag
logFlag = ' ! Set by {}\n'.format(thisFileName)

# Small functions that are only useful here
def convCheck(dataFile, kill=True):
    try:
        f = checkConvergence(dataFile)
        return f
    except (IOError, KeyError, ValueError):
        if kill:
            errStr = 'It appears that the calculation that produced {} did not converge correctly. '.format(dataFile)
            errStr += 'Only properly-converged calculations can be used to spawn Phi1 calculations. '
            errStr += 'Please correct or delete the problematic calculation before trying again. '
            errStr += 'You may wish to check the convergence of all your calculations before running this script.'
            raise IOError(errStr)
        else:
            return None

def makeNewLine(var):
    newLine = '\t' + var + ' = ' + newRunParams[var]['val'] + logFlag
    return newLine

# Loop through the files and process them
for (inDir, outDir) in zip(inDirs, outDirs):

    dataFiles = findFiles('sfincsOutput.h5', inDir, raiseError=True) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in
    dataFileSubdirs = [dirname(dataFile) for dataFile in dataFiles]

    # First, work out which directories have files that need to be copied over
    needToCopy = []
    skip = []
    for i, dataFile in enumerate(dataFiles):

        inSubDir = dataFileSubdirs[i]
        dataDepth = len(dataFile.replace(inDir+'/','').split('/')) # 2 if only radial directories are present, 3 if radial and Er directories are present

        if inSubDir in skip:
            continue

        if dataDepth == 2:
            f = convCheck(dataFile, kill=True)
            needToCopy.append((inSubDir, f))
        
        elif dataDepth == 3: # Must choose proper Er subdirectory to use
            
            radialSubDir = dirname(inSubDir)
            radialMatching = []
            for localDataFile in dataFiles:
                if radialSubDir in localDataFile:
                    f = convCheck(localDataFile, kill=False)
                    if f is not None:
                        radialMatching.append((dirname(localDataFile), f))

            if len(radialMatching) == 0:
                errStr = 'All the calculations that live in the subdirectories of {} appear to have failed. '.format(radialSubDir)
                errStr += 'There must be at least one successful calculation in each radial directory for this script to work.'
                raise IOError(errStr)
            
            extCurs = []
            for matching in radialMatching:
                f = matching[1]
                normalizedAreaFactor = f['VPrimeHat'][()] # = dVHat/dpsiHat
                radialCurrent_vm_psiHat = np.dot(f['Zs'][()], f['particleFlux_vm_psiHat'][()])
                extensiveRadialCurrent = np.abs(normalizedAreaFactor * radialCurrent_vm_psiHat[0])
                extCurs.append(extensiveRadialCurrent)

            minJrInd = np.argmin(extCurs)
            for item in radialMatching:
                skip.append(item[0]) # We don't want to loop into the min(|Jr|) directory in the future
            ErSubDirInfoToUse = radialMatching[minJrInd]
            needToCopy.append(ErSubDirInfoToUse)

        else:
            raise IOError('The structure of the directory {} seems to be irregular.'.format(inDir))
        
    # Now actually copy over files and edit them as needed
    for copyTuple in needToCopy:
        
        copyDir = copyTuple[0]
        sfincsData = copyTuple[1]
        outSubDir = copyDir.replace(inDir, outDir)
        
        # Make target directory if it does not exist
        _ = makeDir(outSubDir) # Note that this script has file overwrite powers!

        # Load in input.namelist file from inDir
        try:
            with open(join(copyDir, 'input.namelist'),'r') as f:
                originalInputLines = f.readlines()
        except FileNotFoundError:
            raise IOError('<sfincsDir> directory(ies) must always contain subdirectory(ies) with "input.namelist" files that were used as SFINCS inputs.')
        
        # Sort out the new entries for input.namelist
        newInputLines = []
        for line in originalInputLines:
            
            # Identify relevant pieces of text from the line
            splitLine = line.split('=')
            var = splitLine[0].strip().split('!ss')[-1].strip()
            if len(splitLine) > 1:
                val = splitLine[1].strip().split('!')[0].strip()
            else:
                val = None
            
            # We'll need the radial coordinate choices for createing a runspec.dat file # FIXME I think you'll still need this for setting Er properly
            if var == 'inputRadialCoordinate':
                inputRadialCoordinate = val
            elif var == 'inputRadialCoordinateForGradients':
                inputRadialCoordinateForGradients = val

            # Fix the line if necessary
            if var in newRunParams.keys():
                newLine = makeNewLine(var)
                newRunParams[var]['used'] = True
            else:
                newLine = line

            # Record the line so it can be written later
            newInputLines.append(newLine)

        # Check if any new settings need to be added
        for key,val in newRunParams.items():
            if val['used'] == False:
                print('FOUND FALSE', key)
                newLine = makeNewLine(key)
                ind = newInputLines.index('&' + val['paramList'] + '\n') + 1
                newInputLines.insert(ind, newLine)
            newRunParams[key]['used'] = False # Must reset before the next file

        # Write new input.namelist file
        stringToWrite = ''.join(newInputLines)
        outFile = join(outSubDir, 'input.namelist')
        writeFile(outFile, stringToWrite)

    # FIXME you need to copy over jobs files too!
    # FIXME perhaps make a 'queue' of files to send to slurm and only do it once all the convergence checks have been passed?

# FIXME don't forget to add in any other options that were missed in the original file!
