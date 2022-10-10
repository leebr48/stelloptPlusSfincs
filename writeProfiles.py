# This script creates a SFINCS-readable profiles file.

# Import necessary packages
import numpy as np
from IO import getArgs, getFileInfo, cleanStrings, listifyBEAMS3DFile, extractDataList, makeProfileNames, generatePreamble, generateDataText
from dataProc import findMinMax, scaleData, nonlinearInterp

# Get command line arguments
args = getArgs()

# Name input and output files
inFile, inFileName, inFilePath, outFilePath = getFileInfo(args.inFile[0], args.saveLoc[0])

outFileName = 'profiles' # Name mandated by SFINCS

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

radii = list(np.linspace(start=radialBounds['min'], stop=radialBounds['max'], num=args.numInterpSurf[0], endpoint=True))

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
