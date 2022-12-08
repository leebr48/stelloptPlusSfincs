# This script is used to copy input.namelist files from radial or electric field subdirectories for calculations that did not include Phi1 and convert them to running with Phi1.
# job.sfincsScan files are also copied.
# If multiple electric field subdirectories are available, this script will pick the one with the lowest |Jr| and only copy those files.

from os.path import dirname, abspath, join
from inspect import getfile, currentframe
from shutil import copy
from subprocess import run
import sys
import numpy as np
thisFile = abspath(getfile(currentframe()))
thisDir = dirname(thisFile)
sys.path.append(join(thisDir, 'src/'))
from IO import getPhi1SetupArgs, getFileInfo, adjustInputLengths, makeDir, findFiles, radialVarDict, writeFile, messagePrinter, saveTimeStampFile
from dataProc import checkConvergence, convertRadDer
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

# Collect some variables for later
newRunParams = {'ambipolarSolve': {'val': '.false.', 'paramList': 'general', 'used': False},
                'includePhi1': {'val': '.true.', 'paramList': 'physicsParameters', 'used': False}
                }

radialVars = radialVarDict()
jobFileName = 'job.sfincsScan'

logFlag = ' ! Set by {}\n'.format(thisFileName)
logFileString = 'The following automation tasks were carried out:\n'

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

def splitLine(line):
    splitLine = line.split('=')
    var = splitLine[0].strip().split('!ss')[-1].strip()
    if len(splitLine) > 1:
        val = splitLine[1].strip().split('!')[0].strip()
    else:
        val = None
    return var, val

def makeNewLine(var):
    newLine = '\t' + var + ' = ' + newRunParams[var]['val'] + logFlag
    return newLine

def ErDefs():
    l = ['dPhiHatd{}'.format(var) for var in list(radialVars.values())[0:-1]]
    l = ['dPhiHatd{}'.format(radialVars[i]) for i in range(4)]
    l.append('Er')
    return l

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

        # Ensure that the loaded run didn't already include a Phi1 run
        trueInt = f['integerToRepresentTrue'][()]
        wasPhi1Used = f['includePhi1'][()]
        if wasPhi1Used == trueInt:
            raise IOError('The run that created {} already included Phi1!'.format(dataFile))
        
    # Now actually copy over files and edit them as needed
    outSubDirs = []
    for copyTuple in needToCopy:
        
        copyDir = copyTuple[0]
        sfincsData = copyTuple[1]
        outSubDir = copyDir.replace(inDir, outDir)
        outSubDirs.append(outSubDir)
        
        # Load in input.namelist file from inDir
        try:
            with open(join(copyDir, 'input.namelist'),'r') as f:
                originalInputLines = f.readlines()
        except FileNotFoundError:
            raise IOError('<sfincsDir> directory(ies) must always contain subdirectory(ies) with "input.namelist" files that were used as SFINCS inputs.')
        
        # Make target directory if it does not exist
        _ = makeDir(outSubDir) # Note that this script has file overwrite powers!
        
        # Sort out the new entries for input.namelist
        newInputLines = []
        for line in originalInputLines:
            
            # Identify relevant pieces of text from the line
            var, val = splitLine(line)
            
            # We'll need the radial coordinate choices for setting Er properly
            if var == 'inputRadialCoordinateForGradients':
                inputRadialCoordinateForGradients = int(val)
            
            # If the electric field is set in the input file, remove that and add it back later
            if var in ErDefs():
                continue

            # Fix the line if necessary
            if var in newRunParams.keys():
                newLine = makeNewLine(var)
                newRunParams[var]['used'] = True
            else:
                newLine = line

            # Record the line so it can be written later
            newInputLines.append(newLine)

        # Set the electric field appropriately
        Er = sfincsData['Er'][()] # Note that if ambipolarSolve was used in the previous run, only Er will have the value determined by the root-finding algorithm.
        # The dPhiHatd* variables will only have their seed values.
        ErID = ErDefs().index('Er')
        aHat = sfincsData['aHat'][()]
        psiAHat = sfincsData['psiAHat'][()]
        psiN = sfincsData['psiN'][()]
        generalizedErVal = convertRadDer(ErID, Er, inputRadialCoordinateForGradients, aHat, psiAHat, psiN, XisPhi=True)
        newRunParams[ErDefs()[inputRadialCoordinateForGradients]] = {'val': str(generalizedErVal), 'paramList': 'physicsParameters', 'used': False}

        # Check if any new settings need to be added
        for key,val in newRunParams.items():
            if val['used'] == False:
                newLine = makeNewLine(key)
                ind = newInputLines.index('&' + val['paramList'] + '\n') + 1
                newInputLines.insert(ind, newLine)
            newRunParams[key]['used'] = False # Must reset before the next file

        # Write new input.namelist file
        stringToWrite = ''.join(newInputLines)
        outFile = join(outSubDir, 'input.namelist')
        writeFile(outFile, stringToWrite, silent=True)

        # Copy the job.sfincsScan file to the new directory
        copy(join(copyDir, jobFileName), outSubDir)

    messagePrinter('All relevant files have been copied from {} to {}.'.format(inDir, outDir))
    logFileString += '\tinput.namelist file(s) were pulled from {}, converted to include Phi1, and written\n'.format(inDir)

    # Now that all the files have been written, send them to Slurm if necessary
    if not args.noRun:

        cmd = ['sbatch', jobFileName]
        for outSubDir in outSubDirs:
            run(cmd, cwd=outSubDir)
    
        messagePrinter('All runs have been submitted for {}.'.format(outDir))
        logFileString += '\tjob(s) were submitted to be run\n'

    # Write a log file
    logFileString += 'at this time:\n\t'
    saveTimeStampFile(outDir, 'automatedSetupLog', logFileString)
