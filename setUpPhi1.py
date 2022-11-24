# FIXME This script is used to copy input.namelist files and convert them to running with Phi1.

from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
thisFile = abspath(getfile(currentframe()))
thisDir = dirname(thisFile)
sys.path.append(join(thisDir, 'src/'))
from IO import getPhi1SetupArgs, getFileInfo, adjustInputLengths
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
newRunParams = {'scanType': {'preFlag': '!ss ', 'val': '21'},
                'ambipolarSolve': {'preFlag': '\t', 'val': '.false.'},
                'includePhi1': {'preFlag': '\t', 'val': '.true.'}
                } # If variable is for sfincsScan, set preFlag to '!ss ', otherwise set it to '\t'

# I/O flag
logFlag = ' ! Set by {}'.format(thisFileName)

# Loop through the files and process them
for inDir in inDirs: 
    
    # Load in input.namelist file from inDir
    try:
        with open(join(inDir, 'input.namelist'),'r') as f:
            originalInputLines = f.readlines()
    except FileNotFoundError:
        raise IOError('<sfincsDir> directory(ies) must always contain an "input.namelist" file that was used in conjunction with sfincsScan.')
    
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
        
        # We'll need the radial coordinate choices for createing a runspec.dat file
        if var == 'inputRadialCoordinate':
            inputRadialCoordinate = val
        elif var == 'inputRadialCoordinateForGradients':
            inputRadialCoordinateForGradients = val

        # Fix the line if necessary
        if var in newRunParams.keys():
            newLine = newRunParams[var]['preFlag'] + var + ' = ' + newRunParams[var]['val'] + logFlag
        else:
            newLine = line

        # Record the line so it can be written later
        newInputLines.append(newLine)

    # FIXME should probably load the h5 files and such in here, maybe even prior to the namelist... and maybe write the namelist after so you can make final touches using data

# FIXME you need to add a dummy electric field value (with the right radial variable) and you need to turn on includePhi1... could you just make all your files be written with these variables for ease? (I think you've done this, check)

# FIXME don't forget to add in any other options that were missed in the original file! -- might only do this if sfincsScan throws an error... If that happens, can perhaps write them for all the files
