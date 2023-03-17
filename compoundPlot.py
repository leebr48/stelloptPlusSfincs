# This script can be used to plot *.dat information created by plot.py. Only one plot can be created at a time.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getCompoundPlotArgs, getFileInfo, makeDir

# Get user inputs
args = getCompoundPlotArgs()

# Small function useful only in this script
def setLength(inList, desiredLength):
    
    if inList == ['']:
        outList = [None] * desiredLength
    elif inList != [''] and len(inList) == 1:
        outList = inList * desiredLength
    elif len(inList) == desiredLength:
        outList = inList
    else:
        raise IOError('If you manually specify <colors>, <lineStyles>, <markers>, or <zorder> you must input either one argument (for all curves) or one argument for each curve. In this case, there are {} curves.'.format(desiredLength))

    return outList

# Regularize some input data
_, _, _, outFilePath, outFile = getFileInfo('arbitrary', args.saveLoc, args.fileName+'.'+args.fileType)

_ = makeDir(outFilePath) # Note that this script has file overwrite powers!

regData = []
for item in args.data:
    temp, _, _, _, _  = getFileInfo(item, '/arbitrary/path/', 'arbitrary')
    regData.append(temp)

# Load input data
IVs = []
DVs = []
lineCount = 0
for item in regData:
    temp = np.loadtxt(item)
    lineCount += temp.shape[1] - 1
    IVs.append(temp[:, 0])
    DVs.append(temp[:, 1:])

# Set visual parameters and perform some input checks
if args.legend != [''] and len(args.legend) != lineCount:
    raise IOError('If you manually specify <legend>, you must input one entry for each curve. In this case, there are {} curves.'.format(lineCount))

colors = setLength(args.colors, lineCount)
lineStyles = setLength(args.lineStyles, lineCount)
markers = setLength(args.markers, lineCount)
zorders = setLength(args.zorders, lineCount)

# Plot
plt.rcParams['font.size'] = str(args.fontSize)

fig, ax = plt.subplots()

styleInd = 0
for IV, subDVs in zip(IVs, DVs):
    
    for DV in subDVs.T:
        
        plotArgs = [IV, DV]
        plotKwargs = {'color':colors[styleInd], 'linestyle':lineStyles[styleInd], 'marker':markers[styleInd], 'zorder':zorders[styleInd]}

        if args.plotType[0] == 'linear':
            ax.plot(*plotArgs, **plotKwargs)
        elif args.plotType[0] == 'semilogy':
            ax.semilogy(*plotArgs, **plotKwargs)

        styleInd += 1

if args.hlines != ['']:
    for yval in args.hlines:
        ax.axhline(y=yval, color='black', linestyle='-', zorder=0)

if args.vlines != ['']:
    for xval in args.vlines:
        ax.axvline(x=xval, color='black', linestyle='-', zorder=0)

ax.set_xlabel(args.xlabel[0])
ax.set_ylabel(args.ylabel[0])

if args.legend != ['']:
    ax.legend(args.legend, loc='best')

if args.xmin[0] != None:
    plt.xlim(xmin=args.xmin[0])

plt.margins(x=0.01)

if args.ymargin != None:
    plt.margins(y=args.ymargin)

ax.xaxis.set_major_locator(MultipleLocator(args.xtick[0]))

plt.tight_layout()
fig.savefig(outFile, bbox_inches='tight', pad_inches=0, dpi=400)
