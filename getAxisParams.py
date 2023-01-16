# This is a one-off script that pulls the axis parameters from a VMEC wout file and prints them in such a way that they can be easily copy-pasted into a VMEC input file.

# Load packages
from scipy.io import netcdf
import argparse
from os.path import isdir
import numpy as np

# Load argument
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--wout', type=str, nargs=1, required=True, help='wout file from which to pull axis information. Note that stellarator symmetry is assumed!')
args = parser.parse_args()

if isdir(args.wout[0]):
    raise IOError('The input given in <wout> must be a file.')

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
