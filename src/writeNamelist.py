# This script creates a SFINCS-readable input.namelist file.

def run(profilesInUse, saveLocUse, eqInUse, bcSymUse):
    
    '''
    The inputs are set by a wrapper script.
    '''

    # Import necessary modules
    from IO import getRunArgs, getFileInfo, cleanStrings, listifyBEAMS3DFile, extractScalarData, radialVarDict, writeFile
    from dataProc import scaleInputData, findNumCalcs

    # Get command line arguments
    args = getRunArgs()

    # Name input and output files
    profilesFile, _, _, _, outFile = getFileInfo(profilesInUse, saveLocUse, 'input.namelist') # Name mandated by SFINCS
    eqFile, _, _, _, _ = getFileInfo(eqInUse, '/arbitrary/path/', 'arbitrary')

    # List out some hard-coded variables
    profilesScheme = 1 # The profile information is specified on many flux surfaces rather than using polynomials, simply because it's easier and we don't need to worry about fit quality as much
    ambipolarSolveOption = 3 # Use a Newton method
    VMECRadialOption = 0 # Interpolate when the target surface does not exactly match a VMEC flux surface
    Delta = str(4.5694e-3).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    alpha = str(1.0e+0).lower().replace('e','d') # Default -> makes reference quantities sensible/easy
    nu_n = str(-1) # Not default... negative value initiates auto-calculation based on conditions of first species (assumed to be electrons)
    collisionOperator = 0 # (Default) Full linearized Fokker-Planck operator
    includeXDotTerm = '.true.' # (Default) Necessary to calculate full trajectories
    includeElectricFieldTermInXiDot = '.true.' # (Default) Necessary to calculate full trajectories
    export_full_f = '.false.' # Whether or not to save the full distribution function in the output file 
    export_delta_f = '.false.' # Whether or not to save the departure from the Maxwellian distribution function in the output file

    # Load necessary variables from profilesFile
    varsOfInterest = cleanStrings(['NI_AUX_M', 'NI_AUX_Z']) # STELLOPT has the mass and charge of electrons built in, so only the ions need to be specified
    listifiedInFile = listifyBEAMS3DFile(profilesFile)
    dataOfInterest = extractScalarData(listifiedInFile, varsOfInterest)
    dataOfInterest['m'].insert(0, args.assumedSpeciesMass[0])
    dataOfInterest['z'].insert(0, args.assumedSpeciesCharge[0])
    scaledData = scaleInputData(dataOfInterest, profiles=False)
    mHatsList = scaledData['m']
    ZsList = scaledData['z']
    
    if len(mHatsList) != len(ZsList):
        raise IOError('It appears that some data is specified incorrectly in {}. The number of species (according to the "Z" and "M" namelist items) must be consistent.'.format(profilesFile))
    numSpecies = len(ZsList)
    
    mHats = ' '.join(['{:.15e}'.format(mHat).replace('e','d') for mHat in mHatsList])
    Zs = ' '.join([str(Z) for Z in ZsList])

    # Sort out some variables set by command line inputs
    if args.resScan:
        scanType = 1
    elif args.numErSubscan[0] == 0:
        scanType = 4
    else:
        scanType = 5

    magneticDriftScheme = args.driftScheme[0] # Whether or not to include angular drifts, and if so, what model to use

    if args.includePhi1:
        includePhi1 = '.true.' # Calculate angular variation of the electric potential on a flux surface
    else:
        includePhi1 = '.false.' # Assume the electric potential is a flux function

    if args.ambiSolve:
        ambipolarSolve = '.true.' # Determine the ambipolar Er
    else:
        ambipolarSolve = '.false.' # Use the given seed value of Er

    Er_search_tolerance_f = str(args.maxRootJr[0]).lower().replace('e','d') # Root-finding tolerance (radial current in SFINCS internal units) - lower than default to help ensure ambipolar fluxes
    
    Er_min = args.minSolverEr[0]
    Er_max = args.maxSolverEr[0]
    
    eqFileExt = eqFile.split('.')[-1]
    if eqFileExt == 'bc' and bcSymUse == 'sym':
        geometryScheme = 11
    elif eqFileExt == 'bc' and bcSymUse == 'asym':
        geometryScheme = 12
    else:
        geometryScheme = 5
        
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
    stringToWrite += '!ss profilesScheme = {} ! How the profile information is specified\n'.format(profilesScheme)
    stringToWrite += '!ss Nradius = {} ! Number of flux surfaces on which to perform full SFINCS calculations if sfincsScan is called appropriately\n'.format(args.numCalcSurf[0])
    stringToWrite += '!ss {}_min = {} ! Lower bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.minRad[0])
    stringToWrite += '!ss {}_max = {} ! Upper bound for the radial scan\n'.format(radialVars[args.radialVar[0]], args.maxRad[0])
    stringToWrite += '\n'

    stringToWrite += '&general\n'
    stringToWrite += '\tambipolarSolve = {} ! Whether or not to determine the ambipolar Er\n'.format(ambipolarSolve)
    stringToWrite += '\tambipolarSolveOption = {} ! Specifies the root-finding algorithm to use\n'.format(ambipolarSolveOption)
    stringToWrite += '\tEr_search_tolerance_f = {} ! Root-finding tolerance (radial current in SFINCS internal units)\n'.format(Er_search_tolerance_f)
    stringToWrite += '\tEr_min = {} ! Minimum value of Er (= -dPhiHatdrHat) accessible to ambipolarSolve.\n'.format(Er_min)
    stringToWrite += '\tEr_max = {} ! Maximum value of Er (= -dPhiHatdrHat) accessible to ambipolarSolve.\n'.format(Er_max)
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&geometryParameters\n'
    stringToWrite += '\tgeometryScheme = {} ! Set how the magnetic geometry is specified\n'.format(geometryScheme)
    stringToWrite += '\tinputRadialCoordinate = {} ! {}\n'.format(args.radialVar[0], selectedRadialVar)
    stringToWrite += '\t{}_wish = {} ! Surface on which to perform the resolution scan (will be overwritten for other applications)\n'.format(selectedRadialVar, args.minRad[0])
    stringToWrite += '\tinputRadialCoordinateForGradients = {} ! {}\n'.format(args.radialGradientVar[0], selectedRadialGradientVar)
    stringToWrite += '\tVMECRadialOption = {} ! Interpolate when the target surface does not exactly match a VMEC flux surface\n'.format(VMECRadialOption)
    stringToWrite += '\tequilibriumFile = "{}"\n'.format(eqFile)
    stringToWrite += '\tmin_Bmn_to_load = {} ! Only Fourier modes of at least this size will be loaded from the equilibriumFile\n'.format(args.minBmn[0])
    if geometryScheme == 5:
        stringToWrite += '\tVMEC_Nyquist_option = {} ! If 2, include the larger poloidal and toroidal mode numbers in the xm_nyq and xn_nyq arrays (where available)\n'.format(args.Nyquist[0])
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&speciesParameters\n'
    stringToWrite += '\tZs = {} ! Charge of each species in units of the proton charge\n'.format(Zs)
    stringToWrite += '\tmHats = {} ! Mass of each species in units of the proton mass\n'.format(mHats)
    stringToWrite += '\tnHats = {} ! Density of each species to use for the resolution scan (may be ignored for other applications)\n'.format(nHats)
    stringToWrite += '\tTHats = {} ! Temperature of each species to use for the resolution scan (may be ignored for other applications)\n'.format(THats)
    stringToWrite += '\tdNHatd{}s = {} ! Radial derivative of density for each species to use for the resolution scan (may be ignored for other applications)\n'.format(selectedRadialGradientVar, dNHatDer)
    stringToWrite += '\tdTHatd{}s = {} ! Radial derivative of temperature for each species to use for the resolution scan (may be ignored for other applications)\n'.format(selectedRadialGradientVar, dTHatDer)
    stringToWrite += '/\n'
    stringToWrite += '\n'

    stringToWrite += '&physicsParameters\n'
    stringToWrite += '\tDelta = {} ! Sets reference units\n'.format(Delta)
    stringToWrite += '\talpha = {} ! Sets reference units\n'.format(alpha)
    stringToWrite += '\tnu_n = {} ! Sets reference units\n'.format(nu_n)
    stringToWrite += '\tcollisionOperator = {} ! Specifies collision operator to use\n'.format(collisionOperator)
    stringToWrite += '\tincludeXDotTerm = {} ! This term is necessary to calculate full trajectories\n'.format(includeXDotTerm)
    stringToWrite += '\tincludeElectricFieldTermInXiDot = {} ! This term is necessary to calculate full trajectories\n'.format(includeElectricFieldTermInXiDot)
    stringToWrite += '\tmagneticDriftScheme = {} ! Whether or not to include tangential drifts, and if so, which model to use\n'.format(magneticDriftScheme)
    stringToWrite += '\tincludePhi1 = {} ! Whether or not to include variation of electric potential on the flux surface\n'.format(includePhi1)
    if args.radialGradientVar[0] != 4:
        stringToWrite += '\tdPhiHatd{} = {} ! Seed value of the radial electric field (proxy) for this flux surface (may be ignored)\n'.format(selectedRadialGradientVar, args.seedEr[0])
    else:
        stringToWrite += '\tEr = {} ! Seed value of the radial electric field for this flux surface (may be ignored)\n'.format(args.seedEr[0])
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
    stringToWrite += '\tsolverTolerance = {} ! Tolerance that specifies convergence for the iterative solver\n'.format(solverTol)
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
    stringToWrite += '\texport_delta_f = {} ! Whether or not to save the departure from the Maxwellian distribution function in the output file\n'.format(export_delta_f)
    stringToWrite += '/\n'

    # Write input.namelist file
    writeFile(outFile, stringToWrite)
