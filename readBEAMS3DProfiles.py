# This script reads profile data from the BEAMS3D section of the STELLOPT namelist and converts it to SFINCS-readable profiles. FIXME description when finished

# Import necessary packages
import argparse
import os
import numpy as np
from IO import listifyBEAMS3DFile, extractDataList, makeProfileNames
from dataProc import linearInterp

# Specify and explain command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='Input file name, with path if necessary.')
parser.add_argument('--outFileName', type=str, nargs=1, required=False, default=None, help='Output file suffix (profile.<outFileName>). Defaults to suffix of input file.')
parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profile.<outFileName>. Defaults to <inFile> location.')
parser.add_argument('--minRad', type=float, nargs=1, required=False, default=0.05, help='Minimum value of the generalized radial coordinate for the scan.')
parser.add_argument('--maxRad', type=float, nargs=1, required=False, default=0.95, help='Maximum value of the generalized radial coordinate for the scan.')
parser.add_argument('--numRad', type=float, nargs=1, required=False, default=16, help='Number of radial surfaces on which to perform the scan.')
parser.add_argument('--noEr', action='store_true', required=False, help='Ignore the scan over the radial electric field.')
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

# Extract the data from the BEAMS3D file.
listifiedInFile = listifyBEAMS3DFile(inFile)

prefixesOfInterest = ['POT', 'NE', 'TI', 'NE', 'TE'] # Note that the order must match the column order of the profiles.xxx file! Repeated prefixes are fine. #FIXME how deal with ZEFF and POT? #FIXME this list could/should probably be an input
if args.noEr:
    del prefixesOfInterest[0] # Relies on the fact that Er data is specified before species data in the profiles.xxx file.

varsOfInterest = makeProfileNames(prefixesOfInterest)
dataOfInterest = extractDataList(listifiedInFile, varsOfInterest)

# Interpolate the data in case the radial lists do not all contain the same points.
interpolatedData = linearInterp(dataOfInterest) #FIXME you might need a non-linear spline if you end up having to take derivatives for the potential...

# Gather the components of the profile.xxx file
radial_coordinate_ID = 1 # Corresponds to normalized toroidal flux, which is the VMEC S.

radii = np.linspace(start=args.minRad, stop=args.maxRad, num=args.numRad, endpoint=True)

if args.noEr:
    zeroList = ['0']*len(radii)
    NErs_vec = zeroList
    generalEr_min_vec = zeroList
    generalEr_max_vec = zeroList
else:
    raise AssertionError('FIXME: I cannot handle radial electric fields yet!')

stringToWrite = ''
