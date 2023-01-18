# This is a one-off script that pulls the axis parameters from a VMEC wout file and prints them in such a way that they can be easily copy-pasted into a VMEC input file.

# Load packages
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from scipy.io import netcdf
import numpy as np

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getAxisParamsArgs 

# Load argument
args = getAxisParamsArgs()

# Read wout file
f = netcdf.netcdf_file(args.wout[0], mode='r', mmap=False)

# Get the axis variables
rax = f.variables['raxis_cc'][()]
zax = f.variables['zaxis_cs'][()]

# Print the axis variables
printR = 'RAXIS_CC = ' + np.array2string(rax, separator=' ', precision=14, max_line_width=100000).replace('[','').replace(']','')
printZ = 'ZAXIS_CS = ' + np.array2string(zax, separator=' ', precision=14, max_line_width=100000).replace('[','').replace(']','')

print(printR)
print(printZ)
