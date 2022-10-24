# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, adjustInputLengths, makeDir, findFiles

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

for i,directory in enumerate(IOlists['sfincsDir']):
    
    # Make target directory if it does not exist
    makeDir(saveDefaultTarget[i]) # Note that this script has file overwrite powers!

    # Retrieve the data
    dataFiles = findFiles('sfincsOutput.h5', directory) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in
    
    data = {}
    for file in dataFiles:
        
        fileNameToUse = file
        while True:
            parentDirOfFile = dirname(fileNameToUse) # Walk up the file directory
            if parentDirOfFile == directory:
                break
            NameOfParentDirOfFile = basename(parentDirOfFile)
            dirTypeIndicator = NameOfParentDirOfFile[0].lower()
            if dirTypeIndicator == 'p' or dirTypeIndicator == 'r':
                # FIXME you're in a radial directory
                pass
            if dirTypeIndicator == 'e' or dirTypeIndicator == 'd':
                #FIXME you're in an e-field directory
                pass





            fileNameToUse = parentDirOfFile
