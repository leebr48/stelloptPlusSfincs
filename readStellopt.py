# FIXME this file reads STELLOPT output and plots various quantities from the optimization so that sense can be made of the outputs.

# Import necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from os import environ

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
stelloptDir = environ['STELLOPT_PATH']
sys.path.append(join(stelloptDir, 'pySTEL'))

from IO import getPyStelArgs, getFileInfo, makeDir
from libstell.stellopt import read_stellopt

# Import arguments
args = getPyStelArgs()

# Set up directories
inFileRaw = args.stelloptOut[0]
inFile, _, inDir, _, _ = getFileInfo(inFileRaw, '/arbitrary/path/', 'arbitrary')

if args.saveLoc[0] is None:
    outDir = join(inDir, 'plots')
else:
    _, _, _, outDir, _ = getFileInfo('arbitrary', args.saveLoc[0], 'arbitrary')

_ = makeDir(outDir) # Note that this script has file overwrite powers!

# Load the input file
data = read_stellopt(inFile)

if not data:
    raise IOError('There was an error opening the file {}.'.format(inFile))
