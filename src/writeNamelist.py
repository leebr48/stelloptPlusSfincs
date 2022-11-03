# FIXME your checks for the lengths of M and Z should be done in Python somehow (before you run SFINCS)... unless of course you just trust your profiles, which you may
# This script creates a SFINCS-readable input.namelist file.

def run(profilesInUse, saveLocUse, eqInUse):
    
    '''
    The inputs are set by a wrapper script.
    '''

    # Import necessary modules
    from IO import getRunArgs, getFileInfo, radialVarDict, writeFile
    from dataProc import findNumCalcs

    # Get command line arguments
    args = getRunArgs()

    # Name input and output files
    profilesFile, _, _, _, outFile = getFileInfo(profilesInUse, saveLocUse, 'input.namelist') # Name mandated by SFINCS

    eqFile, _, _, _, _ = getFileInfo(eqInUse, '/arbitrary/path/', 'arbitrary')

    # List out some hard-coded variables
    profilesScheme = 1 # The profile information is specified on many flux surfaces rather than using polynomials, simply because it's easier and we don't need to worry about fit quality as much
    ambipolarSolve = '.true.' # Determine the ambipolar Er
    ambipolarSolveOption = 2 # (Default) Use a Brent method
    geometryScheme = 5 # Read a VMEC wout file to specify the magnetic geometry
    VMECRadialOption = 0 # Interpolate when the target surface does not exactly match a VMEC flux surface
    Delta = str(4.5694e-3).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    alpha = str(1.0e+0).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    nu_n = str(8.330e-3).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    collisionOperator = 0 # (Default) Full linearized Fokker-Planck operator
    includeXDotTerm = '.true.' # (Default) Necessary to calculate full trajectories
    includeElectricFieldTermInXiDot = '.true.' # (Default) Necessary to calculate full trajectories
    magneticDriftScheme = 1 # Whether or not to include angular drifts, and if so, what model to use
    export_full_f = '.true.' # Save the full distribution function in the output file 

    # Load necessary variables from profilesFile
    # FIXME is it worth making a class that would load the BEAMS3D inputs only once and store them? How would you do that? Would it even be useful?
    Zs = ' '.join([str(Z) for Z in args.Zs]) # FIXME load this!
    mHats = ' '.join(['{:.15e}'.format(mHat).replace('e','d') for mHat in args.mHats]) #FIXME load this!
    # FIXME don't forget to add electrons!
    # FIXME could add command line args for a 'default species' to make things general but still play nice with STELLOPT
    # FIXME actually, don't you have to say that you're using electrons in STELLOPT?
    # FIXME will all the quasineutrality stuff explode if you have only one ion species?
    
    if len(Zs) != len(mHats):
        raise IOError('It appears that some data is specified incorrectly in {}. The number of species (according to the "Z" and "M" namelist items) must be consistent.'.format(profilesFile))
    numSpecies = len(Zs)

    # Sort out some variables set by command line inputs
    if args.resScan:
        scanType = 1
    elif args.numManErScan[0] == 0:
        scanType = 4
    else:
        scanType = 5

    ErDiff = args.maxEr[0] - args.minEr[0]
    Er_min = args.minEr[0] - 0.5 * ErDiff # To ensure SFINCS has an appropriate Er bound if scanning Er
    Er_max = args.maxEr[0] + 0.5 * ErDiff # To ensure SFINCS has an appropriate Er bound if scanning Er
    
    radialVars = radialVarDict()
    selectedRadialVar = radialVars[args.radialVar[0]]
    selectedRadialGradientVar = radialVars[args.radialGradientVar[0]] # Note that option 4 has special Er behavior (see <help> for details)

    nHats = ' '.join(['{:.15e}'.format(nHat).replace('e','d') for nHat in args.defaultDens*numSpecies])
    THats = ' '.join(['{:.15e}'.format(THat).replace('e','d') for THat in args.defaultTemps*numSpecies])
    dNHatDer = ' '.join(['{:.15e}'.format(dnHat).replace('e','d') for dnHat in args.defaultDensDer*numSpecies])
    dTHatDer = ' '.join(['{:.15e}'.format(dTHat).replace('e','d') for dTHat in args.defaultTempsDer*numSpecies])

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

    stringToWrite += '!ss scanType = {}\n'.format(scanType)
    stringToWrite += '!ss profilesScheme = {} ! The profile information is specified on many flux surfaces rather than using polynomials\n'.format(profilesScheme)
    stringToWrite += '!ss Nradius = {} ! Number of radial surfaces on which to perform full SFINCS calculations\n'.format(args.numCalcSurf[0])
    stringToWrite += '!ss {}_min = {} ! Lower bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.minRad[0])
    stringToWrite += '!ss {}_max = {} ! Upper bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.maxRad[0])
    stringToWrite += '\n'

    stringToWrite += '&general\n'
    stringToWrite += '\tambipolarSolve = {} ! Determine the ambipolar Er\n'.format(ambipolarSolve)
    stringToWrite += '\tambipolarSolveOption = {} ! Specifies the root-finding algorithm to use\n'.format(ambipolarSolveOption)
    stringToWrite += '\tEr_min = {} ! Minimum value of Er accessible to ambipolarSolve\n'.format(Er_min)
    stringToWrite += '\tEr_max = {} ! Maximum value of Er accessible to ambipolarSolve\n'.format(Er_max)
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&geometryParameters\n'
    stringToWrite += '\tgeometryScheme = {} ! Read a VMEC wout file to specify the magnetic geometry\n'.format(geometryScheme)
    stringToWrite += '\tinputRadialCoordinate = {} ! {}\n'.format(args.radialVar[0], selectedRadialVar)
    if scanType == 1:
        stringToWrite += '\t{}_wish = {} ! Surface on which to perform the resolution scan\n'.format(selectedRadialVar, args.minRad[0])
    stringToWrite += '\tinputRadialCoordinateForGradients = {} ! {}\n'.format(args.radialGradientVar[0], selectedRadialGradientVar)
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
        stringToWrite += '\tdNHatd{}s = {} ! Derivative of density wrt psiN for each species to use for the resolution scan\n'.format(selectedRadialGradientVar, dNHatDer)
        stringToWrite += '\tdTHatd{}s = {} ! Derivative of temperature wrt psiN for each species to use for the resolution scan\n'.format(selectedRadialGradientVar, dTHatDer)
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
    if args.radialGradientVar[0] != 4:
        stringToWrite += '\tdPhiHatd{} = {} ! Value of the radial electric field (proxy) that will be used for all flux surfaces\n'.format(selectedRadialGradientVar, args.seedEr[0])
    else:
        stringToWrite += '\tEr = {} ! Value of the radial electric field that will be used for all flux surfaces\n'.format(args.seedEr[0])
    stringToWrite += '\tmagneticDriftScheme = {} ! Whether or not to include angular drifts, and if so, what model to use\n'.format(magneticDriftScheme)
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
