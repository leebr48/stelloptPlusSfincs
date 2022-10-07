# This script creates a SFINCS-readable profiles file.

# Import necessary packages
import argparse
import os
import numpy as np
from IO import cleanStrings, listifyBEAMS3DFile, extractDataList, makeProfileNames, generatePreamble, generateDataText
from dataProc import findMinMax, scaleData, nonlinearInterp

# Specify command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='File with relevant profiles written as in STELLOPT, with path if necessary. This script currently reads the BEAMS3D section of the STELLOPT namelist file.')
parser.add_argument("--vars", type=str, nargs='*', required=True, help='''Prefixes of variables to be read, normalized, and written. You should enter each prefix in quotes and put spaces between prefixes. The prefix names are not case sensitive. The density and temperature prefixes should come in the format <'N1' 'T1' 'N2' 'T2' ...> where '1' and '2' often indicate species identifiers (such as 'I' or 'E'). Note that you can write duplicate data by repeating entries. For instance, inputting <'NE' 'TI' 'NE' 'TE'> enforces NI=NE. The order in which the species prefixes are specified should match the species order in input.namelist. If you have potential data to input to calculate the radial electric field, 'POT' can be added anywhere in the list. The potential should give -Er when differentiated with respect to the STELLOPT coordinate S, which is psiN in SFINCS.''')
parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profiles. Defaults to <inFile> location.')
parser.add_argument('--numRad', type=int, nargs=1, required=False, default=[1000], help='Number of radial surfaces on which to calculate and write interpolated profile data. This number should be quite large. Default = 1000.')
parser.add_argument('--numErScan', type=int, nargs=1, required=False, default=[5], help='If a radial electric field scan should occur: number of scans to perform. This parameter will be overwritten if Er data is provided. Default = 5.')
parser.add_argument('--minEr', type=float, nargs=1, required=False, default=[-10], help='If a radial electric field scan should occur: minimum value of the generalized Er variable. Note that you may need to change this to get good results. This parameter will be overwritten if Er data is provided. Default = -10.')
parser.add_argument('--maxEr', type=float, nargs=1, required=False, default=[10], help='If a radial electric field scan should occur: maximum value of the generalized Er variable. Note that you may need to change this to get good results. This parameter will be overwritten if Er data is provided. Default = 10.')
parser.add_argument('--constEr', action='store_true', required=False, help='Assume the radial electric field is constant (as in scanType = 4). Er is set in input.namelist in this case.')
parser.add_argument('--phiBar', type=float, nargs=1, required=False, default=[1], help='Reference electrostatic potential in units of kV. Default = 1.')
parser.add_argument('--nBar', type=float, nargs=1, required=False, default=[1e20], help='Reference density in units of m^(-3). Note that Python "E" notation is equivalent to Fortran "D" notation. Default = 1e20.')
parser.add_argument('--TBar', type=float, nargs=1, required=False, default=[1], help='Reference temperature in units of keV. Default = 1.')
args = parser.parse_args()

# Name input and output files
inFile = os.path.abspath(args.inFile[0])
inFileName = inFile.split('/')[-1]
inFilePath = inFile.replace(inFileName,'')

if args.saveLoc == None:
    outFilePath = inFilePath
else:
    outFilePath = os.path.abspath(args.saveLoc[0])

outFileName = 'profiles'

outFile = outFilePath + '/' + outFileName

# Clean input variable names and do some clerical checks
prefixesOfInterest = cleanStrings(args.vars)

# Extract the data from the BEAMS3D input file
listifiedInFile = listifyBEAMS3DFile(inFile)

varsOfInterest = makeProfileNames(prefixesOfInterest)
dataOfInterest = extractDataList(listifiedInFile, varsOfInterest)

ErDataAvailable = False
if 'pot' in dataOfInterest.keys():
    ErDataAvailable = True

radialBounds = findMinMax(dataOfInterest)

# Scale the data according to the reference variable values.
scaledData = scaleData(dataOfInterest, args.phiBar[0], args.nBar[0], args.TBar[0])

# Interpolate the data in case the radial lists do not all contain the same points.
ders = {}
for key,val in scaledData.items():
    ders[key] = 0

if not args.constEr:
    ders['pot'] = 1 # Only take a derivative when we'll need it for further calculations

interpolatedData = nonlinearInterp(scaledData, ders)

# Gather the components of profiles file
radial_coordinate_ID = 1 # Corresponds to normalized toroidal flux, which is S in STELLOPT and psiN in SFINCS.

radii = list(np.linspace(start=radialBounds['min'], stop=radialBounds['max'], num=args.numRad[0], endpoint=True))

if args.constEr:
    # Note that these quantities must be specified for scanType = 5, but they are ignored for scanType = 4.
    NErs = lambda x: 0
    generalEr_min = lambda x: 0
    generalEr_max = lambda x: 0

else:
    
    if ErDataAvailable:
        # This is scanType = 5, but with a single value of the electric field specified (no scan).
        NErs = lambda x: 1
        generalEr_min = interpolatedData['pot']
        generalEr_max = interpolatedData['pot']

    else:
        # This is scanType = 5 with a proper electric field scan.
        NErs = lambda x: args.numErScan[0]
        generalEr_min = lambda x: args.minEr[0]
        generalEr_max = lambda x: args.maxEr[0]

funcs = [NErs, generalEr_min, generalEr_max]
funcs.extend([interpolatedData[prefix] for prefix in prefixesOfInterest if prefix != 'pot'])

# Get the string to write in profiles file
stringToWrite = generatePreamble(radial_coordinate_ID)
stringToWrite += generateDataText(radii, *funcs)

# Write profiles file
with open(outFile, 'w') as f:
    f.write(stringToWrite)
