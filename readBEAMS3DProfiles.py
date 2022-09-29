# This script reads profile data from the BEAMS3D section of the STELLOPT namelist and converts it to SFINCS-readable profiles. FIXME description when finished

# Import necessary packages
import argparse
from IO import listifyBEAMS3DFile

# Specify and explain command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFile', type=str, nargs=1, required=True, help='Input file name.')
parser.add_argument('--outFile', type=str, nargs=1, required=False, default=None, help='Output file suffix (profile.<out>). Defaults to suffix of input file.') # Output file suffix (profile.<out>)
args = parser.parse_args()

# Name input and output files appropriately in the code
inFile = args.inFile[0].strip()
inFileName = inFile.split('/')[-1]
filesSuffix = inFileName.strip('input.')
inFilePath = inFile.strip(inFileName)

if args.outFile == None:
    outFile = 'profile.'+filesSuffix
else:
    outFile = args.outFile[0]

out = listifyBEAMS3DFile(inFile)
