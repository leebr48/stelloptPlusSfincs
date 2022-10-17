# This script creates a SFINCS-readable input.namelist file.

def run():

    # Import necessary modules
    from IO import getArgs, getFileInfo, writeFile
    from dataProc import findNumCalcs

    # Get command line arguments
    args = getArgs()

    # Name input and output files
    _, _, _, _, outFile = getFileInfo(args.profilesIn[0], args.saveLoc[0], 'input.namelist') # Name mandated by SFINCS

    eqFile, _, _, _, _ = getFileInfo(args.eqIn[0], '/arbitrary/path/', 'arbitrary')

    # List out some hard-coded variables
    profilesScheme = 1 # The profile information is specified on many flux surfaces rather than using polynomials, simply because it's easier and we don't need to worry about fit quality as much
    geometryScheme = 5 # Read a VMEC wout file to specify the magnetic geometry
    inputRadialCoordinateForGradients = 1 # Derivatives wrt psiN (the STELLOPT "S") - unlikely to change because S is used as the independent variable for specifying profiles in STELLOPT
    VMECRadialOption = 0 # Interpolate when the target surface does not exactly match a VMEC flux surface
    Delta = str(4.5694e-3).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    alpha = str(1.0e+0).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    nu_n = str(8.330e-3).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    collisionOperator = 0 # (Default) Full linearized Fokker-Planck operator
    includeXDotTerm = '.true.' # (Default) Necessary to calculate full trajectories
    includeElectricFieldTermInXiDot = '.true.' # (Default) Necessary to calculate full trajectories
    export_full_f = '.true.' # Save the full distribution function in the output file 

    # Sort out some variables prior to string creation
    if args.resScan:
        scanType = 1
    elif args.constEr[0] is not False:
        scanType = 4
    else:
        scanType = 5

    radialVars = {0:'psiHat', 1:'psiN', 2:'rHat', 3:'rN'}
    selectedRadialVar = radialVars[args.radialVar[0]]

    Zs = ' '.join([str(Z) for Z in args.Zs])
    mHats = ' '.join(['{:.15e}'.format(mHat).replace('e','d') for mHat in args.mHats])
    nHats = ' '.join(['{:.15e}'.format(nHat).replace('e','d') for nHat in args.defaultDens])
    THats = ' '.join(['{:.15e}'.format(THat).replace('e','d') for THat in args.defaultTemps])
    dNHatdpsiNs = ' '.join(['{:.15e}'.format(dnHat).replace('e','d') for dnHat in args.defaultDensDer])
    dTHatdpsiNs = ' '.join(['{:.15e}'.format(dTHat).replace('e','d') for dTHat in args.defaultTempsDer])

    NthetaScanVars = findNumCalcs(args.Ntheta[0], args.NthetaScan, powersMode=False)
    NzetaScanVars = findNumCalcs(args.Nzeta[0], args.NzetaScan, powersMode=False)
    NxiScanVars = findNumCalcs(args.Nxi[0], args.NxiScan, powersMode=False)
    NxScanVars = findNumCalcs(args.Nx[0], args.NxScan, powersMode=False)
    NLScanVars = findNumCalcs(args.NL[0], args.NLScan, powersMode=False)
    SolverTolScanVars = findNumCalcs(args.solverTol[0], args.solverTolScan, powersMode=True)
    solverTol = str(args.solverTol[0]).lower().replace('e','d')

    # Create the string to be written
    stringToWrite = '! Input file for SFINCS version 3\n'
    stringToWrite += '\n'

    stringToWrite += '!ss scanType = {} ! Scans of Er are nested within each radial scan\n'.format(scanType)
    stringToWrite += '!ss profilesScheme = {} ! The profile information is specified on many flux surfaces rather than using polynomials\n'.format(profilesScheme)
    stringToWrite += '!ss Nradius = {} ! Number of radial surfaces on which to perform full SFINCS calculations\n'.format(args.numCalcSurf[0])
    stringToWrite += '!ss {}_min = {} ! Lower bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.radialMin[0])
    stringToWrite += '!ss {}_max = {} ! Upper bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.radialMax[0])
    stringToWrite += '\n'

    stringToWrite += '&general\n'
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&geometryParameters\n'
    stringToWrite += '\tgeometryScheme = {} ! Read a VMEC wout file to specify the magnetic geometry\n'.format(geometryScheme)
    stringToWrite += '\tinputRadialCoordinate = {} ! {}\n'.format(args.radialVar[0], selectedRadialVar)
    if scanType == 1:
        stringToWrite += '\t{}_wish = {} ! Surface on which to perform the resolution scan\n'.format(selectedRadialVar, args.radWish[0])
    stringToWrite += '\tinputRadialCoordinateForGradients = {} ! Derivatives wrt psiN (the STELLOPT "S")\n'.format(inputRadialCoordinateForGradients)
    stringToWrite += '\tVMECRadialOption = {} ! Interpolate when the target surface does not exactly match a VMEC flux surface\n'.format(VMECRadialOption)
    stringToWrite += '\tequilibriumFile = "{}"\n'.format(eqFile)
    stringToWrite += '\tmin_Bmn_to_load = {} ! Only Fourier modes of at least this size will be loaded from the equilibriumFile\n'.format(args.minBmn[0])
    stringToWrite += '\tVMEC_Nyquist_option = {} ! If 2, include the larger poloidal and toroidal mode numbers in the xm_nyq and xn_nyq arrays (where available)\n'.format(args.Nyquist[0])
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&speciesParameters\n'
    stringToWrite += '\tZs = {} ! Charge of each species in units of the proton charge\n'.format(Zs)
    stringToWrite += '\tmHats = {} ! Mass of each species in units of the proton mass\n'.format(mHats)
    if scanType == 1:
        stringToWrite += '\tnHats = {} ! Density of each species to use for the resolution scan\n'.format(nHats)
        stringToWrite += '\tTHats = {} ! Temperature of each species to use for the resolution scan\n'.format(THats)
        stringToWrite += '\tdNHatdpsiNs = {} ! Derivative of density wrt psiN for each species to use for the resolution scan\n'.format(dNHatdpsiNs)
        stringToWrite += '\tdTHatdpsiNs = {} ! Derivative of temperature wrt psiN for each species to use for the resolution scan\n'.format(dTHatdpsiNs)
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&physicsParameters\n'
    stringToWrite += '\tDelta = {} ! Default -> makes reference quantities sensible/easy\n'.format(Delta)
    stringToWrite += '\talpha = {} ! Default -> makes reference quantities sensible/easy\n'.format(alpha)
    stringToWrite += '\tnu_n = {} ! Default -> makes reference quantities sensible/easy\n'.format(nu_n)
    stringToWrite += '\tcollisionOperator = {} ! (Default) Full linearized Fokker-Planck operator\n'.format(collisionOperator)
    stringToWrite += '\tincludeXDotTerm = {} ! (Default) Necessary to calculate full trajectories\n'.format(includeXDotTerm)
    stringToWrite += '\tincludeElectricFieldTermInXiDot = {} ! (Default) Necessary to calculate full trajectories\n'.format(includeElectricFieldTermInXiDot)
    # Note that the physics parameters above this point are SFINCS defaults - they are included only for code self-documentation.
    if args.constEr[0] is not False:
        stringToWrite += '\tdPhiHatdpsiN = {} ! Value of the radial electric field (proxy) that will be used for all flux surfaces\n'.format(args.constEr[0])
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&resolutionParameters\n'
    stringToWrite += '\tNtheta = {} ! Number of poloidal grid points (should be odd)\n'.format(args.Ntheta[0])
    stringToWrite += '!ss NthetaMinFactor = {}\n'.format(NthetaScanVars['min'])
    stringToWrite += '!ss NthetaMaxFactor = {}\n'.format(NthetaScanVars['max'])
    stringToWrite += '!ss NthetaNumRuns = {}\n'.format(NthetaScanVars['num'])
    stringToWrite += '\tNzeta = {} ! Number of toroidal grid points per period (should be odd)\n'.format(args.Nzeta[0])
    stringToWrite += '!ss NzetaMinFactor = {}\n'.format(NzetaScanVars['min'])
    stringToWrite += '!ss NzetaMaxFactor = {}\n'.format(NzetaScanVars['max'])
    stringToWrite += '!ss NzetaNumRuns = {}\n'.format(NzetaScanVars['num'])
    stringToWrite += '\tNxi = {} ! Number of Legendre polynomials used to represent the pitch-angle dependence of the distribution function\n'.format(args.Nxi[0])
    stringToWrite += '!ss NxiMinFactor = {}\n'.format(NxiScanVars['min'])
    stringToWrite += '!ss NxiMaxFactor = {}\n'.format(NxiScanVars['max'])
    stringToWrite += '!ss NxiNumRuns = {}\n'.format(NxiScanVars['num'])
    stringToWrite += '\tNx = {} ! Number of grid points in energy used to represent the distribution function\n'.format(args.Nx[0])
    stringToWrite += '!ss NxMinFactor = {}\n'.format(NxScanVars['min'])
    stringToWrite += '!ss NxMaxFactor = {}\n'.format(NxScanVars['max'])
    stringToWrite += '!ss NxNumRuns = {}\n'.format(NxScanVars['num'])
    stringToWrite += '\tNL = {} ! Number of Legendre polynomials used to represent the Rosenbluth potentials\n'.format(args.NL[0])
    stringToWrite += '!ss NLMinFactor = {}\n'.format(NLScanVars['min'])
    stringToWrite += '!ss NLMaxFactor = {}\n'.format(NLScanVars['max'])
    stringToWrite += '!ss NLNumRuns = {}\n'.format(NLScanVars['num'])
    stringToWrite += '\tsolverTolerance = {} ! Tolerance used to define convergence of the iterative solver\n'.format(solverTol)
    stringToWrite += '!ss solverToleranceMinFactor = {}\n'.format(SolverTolScanVars['min'])
    stringToWrite += '!ss solverToleranceMaxFactor = {}\n'.format(SolverTolScanVars['max'])
    stringToWrite += '!ss solverToleranceNumRuns = {}\n'.format(SolverTolScanVars['num'])
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&otherNumericalParameters\n'
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&preconditionerOptions\n'
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&export_f\n'
    stringToWrite += '\texport_full_f = {} ! Whether or not to save the full distribution function in the output file\n'.format(export_full_f)
    stringToWrite += '/\n'

    # Write input.namelist file
    writeFile(outFile, stringToWrite)