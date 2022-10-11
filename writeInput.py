# This script creates a SFINCS-readable input.namelist file.

# Import necessary modules
from IO import getArgs, getFileInfo

# Get command line arguments
args = getArgs()

# Name input and output files
_, _, _, outFilePath = getFileInfo(args.profilesIn[0], args.saveLoc[0])

outFileName = 'input.namelist' # Name mandated by SFINCS

outFile = outFilePath + '/' + outFileName

eqFile, _, _, _ = getFileInfo(args.eqIn[0], '/arbitrary/path/')

# Sort out some variables prior to printing
if args.resScan:
    scanType = 1
elif args.constEr[0]:
    scanType = 4
else:
    scanType = 5

radialVars = {0:'psiHat', 1:'psiN', 2:'rHat', 3:'rN'}

Zs = ' '.join([str(Z) for Z in args.Zs])
mHats = ' '.join(['{:.15e}'.format(mHat).replace('e','d') for mHat in args.mHats])

solverTol = str(args.solverTol[0]).replace('e','d')
    
# Create the string to be printed
#FIXME after getting the standard parameters in here, probably exhaustively set all of them just in case (don't rely on defaults)
stringToWrite = '! Input file for SFINCS version 3\n'
stringToWrite += '\n'

stringToWrite += '!ss scanType = {} ! Scans of Er are nested within each radial scan\n'.format(scanType)
stringToWrite += '!ss profilesScheme = 1 ! The profile information is specified on many flux surfaces rather than using polynomials\n'
stringToWrite += '!ss Nradius = {} ! Number of radial surfaces on which to perform full SFINCS calculations\n'.format(args.numCalcSurf[0])
stringToWrite += '!ss {}_min = {} ! Lower bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.radialMin[0])
stringToWrite += '!ss {}_max = {} ! Upper bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.radialMax[0])
stringToWrite += '\n'

stringToWrite += '&general\n'
#FIXME probably fill out the defaults here
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&geometryParameters\n'
stringToWrite += '\tgeometryScheme = 5 ! Read a VMEC wout file to specify the magnetic geometry\n'
stringToWrite += '\tinputRadialCoordinate = {} ! {}\n'.format(args.radialVar[0], radialVars[args.radialVar[0]])
stringToWrite += '\tinputRadialCoordinateForGradients = 1 ! Derivatives wrt psiN (the STELLOPT "S")\n'
stringToWrite += '\tVMECRadialOption = 0 ! Interpolate when the target surface does not exactly match a VMEC flux surface\n'
stringToWrite += '\tequilibriumFile = "{}"\n'.format(eqFile)
stringToWrite += '\tmin_Bmn_to_load = {} ! Only Fourier modes of at least this size will be loaded from the equilibriumFile\n'
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&speciesParameters\n'
stringToWrite += '\tZs = {} ! Charge of each species in units of the proton charge\n'.format(Zs)
stringToWrite += '\tmHats = {} ! Mass of each species in units of the proton mass\n'.format(mHats)
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&physicsParameters\n'
stringToWrite += '\tDelta = 4.5694d-3 ! Default -> makes reference quantities sensible/easy\n'
stringToWrite += '\talpha = 1d+0 ! Default -> makes reference quantities sensible/easy\n'
stringToWrite += '\tnu_n = 8.330d-3 ! Default -> makes reference quantities sensible/easy\n'
if args.constEr[0]:
    stringToWrite += '\tdPhiHatdpsiN = {} ! Value of the radial electric field (proxy) that will be used for all flux surfaces\n'.format(args.constEr[0])
stringToWrite += '\tcollisionOperator = 0 ! Full linearized Fokker-Planck operator\n'
stringToWrite += '\tincludeXDotTerm = .true. ! Necessary to calculate full trajectories\n'
stringToWrite += '\tincludeElectricFieldTermInXiDot = .true. ! Necessary to calculate full trajectories\n'
stringToWrite += '/\n'
stringToWrite += '\n'

#FIXME maybe don't hard-code scan parameters like this? Maybe use lists in args for brevity? Maybe a standard formula?
stringToWrite += '&resolutionParameters\n'
stringToWrite += '\tNtheta = {} ! Number of poloidal grid points (should be odd)\n'.format(args.Ntheta[0])
stringToWrite += '!ss NthetaMinFactor = 0.5\n'
stringToWrite += '!ss NthetaMaxFactor = 2.5\n'
stringToWrite += '!ss NthetaNumRuns = 35\n'
stringToWrite += '\tNzeta = {} ! Number of toroidal grid points per period (should be odd)\n'.format(args.Nzeta[0])
stringToWrite += '!ss NzetaMinFactor = 0.5\n'
stringToWrite += '!ss NzetaMaxFactor = 2.5\n'
stringToWrite += '!ss NzetaNumRuns = 85\n'
stringToWrite += '\tNxi = {} ! Number of Legendre polynomials used to represent the pitch-angle dependence of the distribution function\n'.format(args.Nxi[0])
stringToWrite += '!ss NxiMinFactor = 0.5\n'
stringToWrite += '!ss NxiMaxFactor = 2.5\n'
stringToWrite += '!ss NxiNumRuns = 175\n'
stringToWrite += '\tNx = {} ! Number of grid points in energy used to represent the distribution function\n'.format(args.Nx[0])
stringToWrite += '!ss NxMinFactor = 0.5\n'
stringToWrite += '!ss NxMaxFactor = 2.5\n'
stringToWrite += '!ss NxNumRuns = 12\n'
stringToWrite += '\tsolverTolerance = {} ! Tolerance used to define convergence of the iterative solver\n'.format(solverTol)
stringToWrite += '!ss solverToleranceMinFactor = 0.1\n'
stringToWrite += '!ss solverToleranceMaxFactor = 10\n'
stringToWrite += '!ss solverToleranceNumRuns = 3\n'
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&otherNumericalParameters\n'
#FIXME add anything?
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&preconditionerOptions\n'
#FIXME add anything?
stringToWrite += '/\n'
stringToWrite += '\n'

stringToWrite += '&export_f\n'
stringToWrite += '\texport_full_f = .true.\n'
stringToWrite += '\texport_delta_f = .true.\n'
stringToWrite += '/\n'

# Write input.namelist file
with open(outFile, 'w') as f:
    f.write(stringToWrite)
