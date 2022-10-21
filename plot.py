# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
from os import makedirs

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, getFileInfo, findFiles

# Get command line arguments
args = getPlotArgs()

# Sort out SFINCS directories
inDirs = [] #FIXME do you need to do this, or can you pass in saveLoc and the input directories and let getFileInfo sort it out?
for unRegDir in args.sfincsDir:
    _, _, _, inDir, _ = getFileInfo('/arbitrary/path/', unRegDir, 'arbitrary')
    inDirs.append(inDir)

for directory in inDirs:
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
