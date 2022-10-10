# This script creates a SFINCS-readable input.namelist file.

# Get command line arguments
args = getArgs()

# Name input and output files
inFile, inFileName, inFilePath, outFilePath = getFileInfo(args.inFile[0], args.saveLoc[0])

outFileName = 'input.namelist' # Name mandated by SFINCS

outFile = outFilePath + '/' + outFileName

if args.constEr:
    scanType = 4
    #FIXME you need to be able to set the electric field value in this case, perhaps manually
else:
    scanType = 5

#FIXME after getting the standard parameters in here, probably exhaustively set all of them just in case (don't rely on defaults)
stringToWrite = '! Input file for SFINCS version 3.\n'
stringToWrite += '\n'
stringToWrite += '!ss scanType = {} ! Scans of Er are nested within each radial scan.\n'.format(str(scanType))
stringToWrite += '!ss profilesScheme = 1 ! The profile information is specified on many flux surfaces rather than using polynomials.\n'
stringToWrite += '!ss Nradius = {} ! Number of radial surfaces on which to perform full SFINCS calculations.\n'.format(str(args.numCalcSurf[0]))
#FIXME you should make the radial variables you use general.
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
stringToWrite +=
