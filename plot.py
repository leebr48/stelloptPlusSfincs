# This script generates plots of interest given a sfincs .h5 file. #FIXME WORK IN PROGRESS! fix all comments and such when finished

# Import necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
from os import makedirs, walk

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getPlotArgs, getFileInfo

# Get command line arguments
args = getPlotArgs()

def findFiles(name, path): #FIXME probably put in IO.py
    result = []
    for root, dirs, files in walk(path):
        if name in files:
            result.append(join(root, name))
    return result

# Sort out SFINCS directories
inDirs = []
for unRegDir in args.sfincsDir:
    _, _, _, inDir, _ = getFileInfo('/arbitrary/path/', unRegDir, 'arbitrary')
    inDirs.append(inDir)
print(inDirs)
quit()
for directory in inDirs:
    dataFiles = findFiles('sfincsOutput.h5', directory) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in
    
    data = {}
    for file in dataFiles:
        
        fileNameToUse = file
        while True:
            parentDirOfFile = dirname(fileNameToUse)
            print(parentDirOfFile)
            print(directory)
            quit()
            if parentDirOfFile == directory:
                break
            NameOfParentDirOfFile = basename(parentDirOfFile)
            fileNameToUse = parentDirOfFile
            print(parentDirOfFile)
