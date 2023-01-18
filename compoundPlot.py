#FIXME this script can be used to plot *.dat information created by plot.py

from os.path import dirname, abspath, join
from inspect import getfile, currentframe

import numpy as np
import argparse 
import pathlib as pl #FIXME kill if you can
import os # FIXME be more specific - probably use entry above
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Read user inputs # FIXME move to IO
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#FIXME you should probably load the Z's as well! or wait... maybe not?
parser.add_argument('--data', nargs='+', required=True, help='*.dat files to be plotted. The first column of each file must be horizontal coordinate values, while the remaining columns must be corresponding vertical coordinate (data) values. The file plot.py produces properly-structured data files automatically. Note that you should choose *.dat files with the same horizontal coordinate. Each column (except the first) in the first file passed to this argument will become a curve in the output plot, then each column (except the first) in the second file, and so on -- this is how one can specify the order of the <legend> and <colors> arguments, for example.')
parser.add_argument('--plotType', type=str, nargs=1, required=False, default=['linear'], help='Type of vertical axis. Options are "linear" and "semilogy".')
parser.add_argument('--xlabel', type=str, nargs=1, required=False, default=[''], help='Label for horizontal axis. Be sure to write in quotes!')
parser.add_argument('--ylabel', type=str, nargs=1, required=False, default=[''], help='Label for vertical axis. Be sure to write in quotes!')
parser.add_argument('--xtick', type=float, nargs=1, required=False, default=[1], help='Sets spacing of the horizontal ticks.')
parser.add_argument('--xmin', type=float, nargs=1, required=False, default=[None], help='Sets minimum value on the horizontal axis.')
parser.add_argument('--legend', nargs='+', required=False, default=[''], help='Legend entries. Can accept one or more arguments. Be sure to write each argument in quotes!')
parser.add_argument('--sourcedir', type=str, required=False, default='') #FIXME use this library's standards
parser.add_argument('--outdir', type=str, required=False, default='') #FIXME use this library's standards
parser.add_argument('--fileName', type=str, required=False, default='composite') #FIXME use this library's standards
parser.add_argument('--fileType', type=str, required=False, default='pdf', help='Type of file you wish to save. Options are "pdf" and "png".')
parser.add_argument('--fontSize', type=float, required=False, default=18, help='Font size for plots.')
parser.add_argument('--ymargin', type=float, required=False, default=None, help='Set vertical axis autoscaling margin. Must be between 0 and 1.')
parser.add_argument('--colors', nargs='+', required=False, default=[''], help='Curve colors. For possible options, see the matplotlib documentation. Note that the default matplotlib colors are from the Tableau palette -- for instance, the default blue is "tab:blue", the default orange is "tab:orange", and so forth. Can accept one or more arguments. Be sure to write each argument in quotes! If one argument is input, all curves will have the same color. If multiple arguments are input, this must be the same as the number of lines that are to be plotted.')
args = parser.parse_args()

if args.plotType[0] not in ['linear', 'semilogy']:
    raise IOError('<plotType> must be either "linear" or "semilogy"')

if args.fileType not in ['pdf', 'png']:
    raise IOError('<fileType> must be either "pdf" or "png"')

if not (0 <= args.ymargin <= 1):
    raise IOError('<ymargin> must be between 0 and 1.')

if not (args.colors == [''] or len(args.colors) == 1 or len(args.colors) == len(args.data)): # FIXME won't work because of how data is stored...
    raise IOError('If you manually specify <colors>, you must input either one color (for all curves) or one color for each curve.')

if args.colors == ['']:
    colors = [None] * len(args.data) # FIXME this won't work because multiple curves are stored in each file...
elif args.colors != [''] and len(args.colors) == 1 and len(args.data) != 1:
    colors = args.colors * len(args.data) # FIXME same problem as above...
else:
    colors = args.colors

# FIXME probably load everything up here and then make IV and DV lists

outdir = str(pl.Path.cwd().joinpath(args.outdir).resolve()) # FIXME use your functions to do this if you can...
sourcedir = str(pl.Path.cwd().joinpath(args.sourcedir).resolve())

if not os.path.exists(outdir):
    os.mkdir(outdir)

plt.rcParams['font.size'] = str(args.fontSize)

# Load variables
loaded = []
for item in args.data:
    temp = np.load(str(pl.Path(sourcedir,item)))
    loaded.append(temp)

# Plot
fig, ax = plt.subplots()

for i, item in enumerate(loaded):
    IV = item[:, 0]
    DV = item[:, 1:]
    
    if args.plotType[0] == 'linear':
        ax.plot(IV, DV, c=colors[i])
    elif args.plotType[0] == 'semilogy':
        ax.semilogy(IV, DV, marker='o', c=colors[i])

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

fig.savefig(str(pl.Path(outdir).joinpath(args.fileName+'.'+args.fileType)), bbox_inches='tight', pad_inches=0, dpi=400)
