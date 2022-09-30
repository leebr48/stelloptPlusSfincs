# This script reads profile data from the BEAMS3D section of the STELLOPT namelist and converts it to SFINCS-readable profiles. FIXME description when finished

# Import necessary packages
import argparse
import os
from IO import listifyBEAMS3DFile, extractDataList, makeProfileNames
from dataProc import linearInterp

# Specify and explain command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='Input file name, with path if necessary.')
parser.add_argument('--outFileName', type=str, nargs=1, required=False, default=None, help='Output file suffix (profile.<outFileName>). Defaults to suffix of input file.')
parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profile.<outFileName>. Defaults to <inFile> location.')
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

prefixesOfInterest = ['NE', 'TE', 'TI', 'ZEFF', 'POT'] #FIXME how deal with ZEFF and POT?
varsOfInterest = makeProfileNames(prefixesOfInterest)
dataOfInterest = extractDataList(listifiedInFile, varsOfInterest)

# FIXME we need to do some linear interpolation... Note that you DO NOT need to change the r coordinate if you use option 1 in sfincdScan (psiN). Note also that you might need a non-linear spline if you end up having to take derivatives for the potential...

# Interpolate the data in case the radial lists do not all contain the same points.
interpolatedData = linearInterp(dataOfInterest)
