# This is a one-off script that pulls the axis parameters from a VMEC wout file and prints them in such a way that they can be easily copy-pasted into a VMEC input file.
# This can be useful for restarting STELLOPT optimizations.

# Load packages
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from scipy.io import netcdf_file
import numpy as np

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getAxisParamsArgs, makeStringForStellopt, messagePrinter

# Load argument
args = getAxisParamsArgs()

# Read wout file
f = netcdf_file(args.wout[0], mode='r', mmap=False)

# Get the axis variables
rax = f.variables['raxis_cc'][()]
zax = f.variables['zaxis_cs'][()]

# Print the axis variables
printR = makeStringForStellopt('RAXIS_CC', rax)
printZ = makeStringForStellopt('ZAXIS_CS', zax)

msg = 'The Fourier coefficients for the magnetic axis are:\n'
msg += printR
msg += printZ
messagePrinter(msg)
