# This script reads profile data from the BEAMS3D section of the STELLOPT namelist and converts it to SFINCS-readable profiles. FIXME description when finished #FIXME all comments and descriptions in this script

# Import necessary packages
import argparse
import os
import numpy as np
from IO import cleanStrings, listifyBEAMS3DFile, extractDataList, makeProfileNames, generatePreamble, generateDataText
from dataProc import scaleData, nonlinearInterp

# Specify and explain command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='Input file name, with path if necessary.')
parser.add_argument("--vars", type=str, nargs='*', required=True, help='''Prefixes of variables to be read, normalized, and written, IN ORDER. You should enter each variable in quotes and put spaces between variables. The input names are not case sensitive. If Er should be included in the calculation, 'POT' should come first. If Er is not included, do not include 'POT'. Whether or not 'POT' is included, the density and temperature data should come in the format <'N1' 'T1' 'N2' 'T2' ...> where '1' and '2' often indicate species identifiers (such as 'I' or 'E'). Note that you can write duplicate data by repeating entries. For instance, inputting <'NE' 'TI' 'NE' 'TE'> enforces NI=NE. The order in which the species information is specified should match that in input.namelist.''')
parser.add_argument('--outFileName', type=str, nargs=1, required=False, default=None, help='Output file suffix (profile.<outFileName>). Defaults to suffix of input file.')
parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profile.<outFileName>. Defaults to <inFile> location.')
parser.add_argument('--minRad', type=float, nargs=1, required=False, default=[0.05], help='Minimum value of the generalized radial coordinate for the scan.')
parser.add_argument('--maxRad', type=float, nargs=1, required=False, default=[0.95], help='Maximum value of the generalized radial coordinate for the scan.')
parser.add_argument('--numRad', type=float, nargs=1, required=False, default=[16], help='Number of radial surfaces on which to perform the scan.')
parser.add_argument('--noEr', action='store_true', required=False, help='Ignore the scan over the radial electric field.')
parser.add_argument('--phiBar', type=float, nargs=1, required=False, default=[1], help='Reference electrostatic potential in units of kV.')
parser.add_argument('--nBar', type=float, nargs=1, required=False, default=[1e20], help='Reference density in units of m^(-3). Note that Python "E" notation is equivalent to Fortran "D" notation.')
parser.add_argument('--TBar', type=float, nargs=1, required=False, default=[1], help='Reference temperature in units of keV.')
args = parser.parse_args()

# Name input and output files
inFile = os.path.abspath(args.inFile[0])
inFileName = inFile.split('/')[-1]
filesSuffix = inFileName.replace('input.','')
inFilePath = inFile.replace(inFileName,'')

if args.saveLoc == None:
    outFilePath = inFilePath
else:
    outFilePath = os.path.abspath(args.saveLoc[0])

if args.outFileName == None:
    outFileName = 'profile.' + filesSuffix
else:
    outFileName = 'profile.' + args.outFileName[0]

outFile = outFilePath + '/' + outFileName

# Get clean input variable names
prefixesOfInterest = cleanStrings(args.vars) #FIXME how deal with ZEFF and POT?

if args.noEr:
    for item in prefixesOfInterest:
        if item.lower() == 'pot':
            raise IOError('If you are not calculating the electric field, you should not specify the potential.')

# Extract the data from the BEAMS3D input file and scale it.
listifiedInFile = listifyBEAMS3DFile(inFile)

varsOfInterest = makeProfileNames(prefixesOfInterest)
dataOfInterest = extractDataList(listifiedInFile, varsOfInterest)

# Scale the data according to the reference variable values.
scaledData = scaleData(dataOfInterest, args.phiBar[0], args.nBar[0], args.TBar[0])

# Interpolate the data in case the radial lists do not all contain the same points.
interpolatedData = nonlinearInterp(scaledData)

# Gather the components of profile.xxx
radial_coordinate_ID = 1 # Corresponds to normalized toroidal flux, which is the VMEC S.

radii = list(np.linspace(start=args.minRad[0], stop=args.maxRad[0], num=args.numRad[0], endpoint=True))

if args.noEr:
    NErs = lambda x: 0
    generalEr_min = lambda x: 0
    generalEr_max = lambda x: 0
    # Note that these quantities only must be specified for scanType = 5. They are ignored if scanType = 4.
    funcs = [NErs, generalEr_min, generalEr_max]
else:
    raise AssertionError('FIXME: I cannot handle radial electric fields yet!')
    #FIXME you need to remake the Er vector stuff here!

otherFuncs = [interpolatedData[prefix] for prefix in prefixesOfInterest]

funcs.extend(otherFuncs)

# Get the string to write in profile.xxx
stringToWrite = generatePreamble(radial_coordinate_ID)
stringToWrite += generateDataText(radii, *funcs)

# Write profile.xxx
with open(outFile, 'w') as f:
    f.write(stringToWrite)
