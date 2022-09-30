# This script reads profile data from the BEAMS3D section of the STELLOPT namelist and converts it to SFINCS-readable profiles. FIXME description when finished

# Import necessary packages
import argparse
import os
from IO import listifyBEAMS3DFile, extractDataList, makeProfileNames

# Specify and explain command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='Input file name, with path if necessary.')
parser.add_argument('--outFileName', type=str, nargs=1, required=False, default=None, help='Output file suffix (profile.<outFileName>). Defaults to suffix of input file.')
parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profile.<outFileName>. Defaults to <inFile> location.')
args = parser.parse_args()

# Name input and output files
#inFile = args.inFile[0].strip()
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
