# This file contains IO helper functions.

def getRunArgs():

    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for writing files and running SFINCS.
    '''

    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--profilesIn', type=str, nargs='*', required=True, help='File(s) with relevant profiles written as in STELLOPT, with path(s) if necessary. This script currently reads the BEAMS3D section of STELLOPT namelist files. Note that you must specify densities rather than Zeff. If you input multiple files, order matters!')
    parser.add_argument('--eqIn', type=str, nargs='*', required=True, help='File(s) from which to load the magnetic equilibrium(ia). Can be VMEC wout file(s) in netCDF or ASCII format, or IPP .bc file(s). If you input multiple files, order matters!')
    parser.add_argument('--bcSymmetry', type=str, nargs='*', required=False, default=['sym'], help='If one or more *.bc files are input via <eqIn>, this setting will control whether SFINCS assumes them to be stellarator-symmetric ("sym") or stellarator-asymmetric ("asym"). If one argument is specified, it will be used for all the <eqIn> files. Note that the length of this argument must be either 1 or equivalent to the length of <eqIn>, even if <eqIn> contains a mix of *.bc and VMEC wout files.')
    parser.add_argument('--minBmn', type=float, nargs=1, required=False, default=[0.0], help='Only Fourier modes of at least this size will be loaded from the <eqIn> file(s).')
    parser.add_argument('--Nyquist', type=int, nargs=1, required=False, default=[2], help='This parameter is only relevant if you are loading VMEC equilibria: include the larger poloidal and toroidal mode numbers in the xm_nyq and xn_nyq arrays, where available, if this parameter is set to 2, and exclude these mode numbers if this parameter is set to 1.')
    parser.add_argument('--numInterpSurf', type=int, nargs=1, required=False, default=[1000], help='Number of radial surfaces on which to calculate and write interpolated profile data. This number should be quite large.')
    parser.add_argument('--radialVar', type=int, nargs=1, required=False, default=[3], help='ID of the radial coordinate used in the input.namelist file to specify which surfaces should be scanned over. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "s"), 2 = rHat, and 3 = rN (which is the STELLOPT "rho")')
    parser.add_argument('--radialGradientVar', type=int, nargs=1, required=False, default=[4], help='ID of the radial coordinate used to take derivatives. Relevant for the generalEr_* parameters in the profiles file and specifying the density and temperature derivatives on a single flux suface. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "s"), 2 = rHat, 3 = rN (which is the STELLOPT "rho"), and 4 = rHat (like option 2, except that Er is used in place of dPhiHatdrHat). The default is recommended.')
    parser.add_argument('--numCalcSurf', type=int, nargs=1, required=False, default=[16], help='Number of radial surfaces on which to perform full SFINCS calculations.')
    parser.add_argument('--minRad', type=float, nargs=1, required=False, default=[0.15], help='Lower bound for the radial scan. If <resScan> is used, the flux surface specified by this parameter will be used for the convergence scan. Note that VMEC has resolution issues near the magnetic axis and SFINCS often converges much slower there due to the relatively low collisionality, so setting <minRad> to be very small may cause problems. If the innermost surface of a loaded equilibrium is outside <minRad>, SFINCS will give nonphysical (usually divergent) answers.')
    parser.add_argument('--maxRad', type=float, nargs=1, required=False, default=[0.95], help='Upper bound for the radial scan.')
    parser.add_argument('--ambiSolve', action='store_true', default=False, help='Enable ambipolarSolve. SFINCS will start from a "seed" value of Er specified using other commands and attempt to modify it such that the radial current is driven to zero (which is what we would expect in a real device). Note that ambipolarSolve searches for *a* root, but it might not find the *correct* root. The script chooseErs.py can help with that.')
    parser.add_argument('--maxRootJr', type=float, nargs=1, required=False, default=[1.0e-12], help='Maximum radial current (defined as in SFINCS) that may be present for a given electric field value to be considered a "root". The default is recommended.')
    parser.add_argument('--loadPot', action='store_true', default=False, help='Load a potential from <profilesIn>. This will overwrite <seedEr>. If you use this option, you must set <numErSubscan> >=1 and <radialGradientVar> = 1. The former requirement ensures the software knows whether or not you wish to use the given potential alone or a range around it, and the latter requirement is required because STELLOPT always specifies the potential profile in terms of "s".')
    parser.add_argument('--seedEr', type=float, nargs=1, required=False, default=[0], help="Input an initial guess for the radial electric field in units of <radialGradientVar>. If you use ambipolarSolve, this seed value may influence whether SFINCS converges to the ion or electron root. This parameter will be overwritten if you trigger an electric field scan with <numErSubscan>.")
    parser.add_argument('--numErSubscan', type=int, nargs=1, required=False, default=[0], help='Number of radial electric field scans to perform within each radial directory. This parameter generates equidistant radial electric field seed values between <minSeedEr> and <maxSeedEr>. This parameter will be overwritten if <resScan> is activated.')
    parser.add_argument('--minSeedEr', type=float, nargs=1, required=False, default=[-5], help='If <loadPot> is used, this value will be added to the values of the loaded potential to determine the minimum seed value of the radial electric field on each flux surface in units of <radialGradientVar>. (Note that for typicaly usage, this value should probably be negative.) If <loadPot> is not used, this parameter gives the mimimum seed value of the radial electric field in units of <radialGradientVar>. You may need to change this parameter to get good results.')
    parser.add_argument('--maxSeedEr', type=float, nargs=1, required=False, default=[5], help='If <loadPot> is used, this value will be added to the values of the loaded potential to determine the maximum seed value of the radial electric field on each flux surface in units of <radialGradientVar>. (Note that for typicaly usage, this value should probably be positive.) If <loadPot> is not used, this parameter gives the maximum seed value of the radial electric field in units of <radialGradientVar>. You may need to change this parameter to get good results.')
    parser.add_argument('--minSolverEr', type=float, nargs=1, required=False, default=[-100], help='Explicitly set the minimum Er (=-dPhiHatdrHat, regardless of <radialGradientVar>) available to ambipolarSolve. This will seldom need to be modified. It is included because the Newton method used by ambipolarSolve can sometimes "get lost" if it is seeded poorly and specify progressively larger |Er| values during the root search. Setting this parameter and <maxSolverEr> closer to the electric field seed value would make the runs fail faster in such situations and therefore save time.')
    parser.add_argument('--maxSolverEr', type=float, nargs=1, required=False, default=[100], help='Explicitly set the maximum Er (=-dPhiHatdrHat, regardless of <radialGradientVar>) available to ambipolarSolve. This will seldom need to be modified. It is included because the Newton method used by ambipolarSolve can sometimes "get lost" if it is seeded poorly and specify progressively larger |Er| values during the root search. Setting this parameter and <minSolverEr> closer to the electric field seed value would make the runs fail faster in such situations and therefore save time.')
    parser.add_argument('--resScan', action='store_true', default=False, help='Triggers a SFINCS resolution scan run.')
    parser.add_argument('--defaultDens', type=float, nargs=1, required=False, default=[1], help='If <resScan> is used, this sets the density of each species in units of 1e20 m^-3. The exact value is probably not important. This setting could break quasineutrality, but again, this is probably not a huge problem for resolution scans. If you are concerned about it, you can edit the input.namelist file directly.')
    parser.add_argument('--defaultTemps', type=float, nargs=1, required=False, default=[1], help='If <resScan> is used, this sets the temperature of each species in keV. The exact value is probably not important.')
    parser.add_argument('--defaultDensDer', type=float, nargs=1, required=False, default=[-0.5e0], help='If <resScan> is used, this sets the derivative of the density of each species (in units of 1e20 m^-3) with respect to the radial variable specified by <radialGradientVar>. The exact value is probably not important.')
    parser.add_argument('--defaultTempsDer', type=float, nargs=1, required=False, default=[-2e0], help='If <resScan> is used, this sets the derivative of the temperature of each species (in keV) with respect to the radial variable specified by <radialGradientVar>. The exact value is probably not important.')
    parser.add_argument('--assumedSpeciesMass', type=float, nargs=1, required=False, default=[9.109383701528e-31], help='Mass in kg of a species that will be included in the SFINCS run(s) regardless of what other species the STELLOPT profiles file specifies. Due to STELLLOPT and SFINCS defaults, this should always be electrons. This option is included for completeness, but you should probably never change it.')
    parser.add_argument('--assumedSpeciesCharge', type=float, nargs=1, required=False, default=[-1.0], help='Charge in units of the proton charge of a species that will be included in the SFINCS run(s) regardless of what other species the STELLOPT profiles file specifies. Due to STELLLOPT and SFINCS defaults, this should always be electrons. This option is included for completeness, but you should probably never change it.')
    parser.add_argument('--Nzeta', type=int, nargs=1, required=False, default=[115], help='Number of toroidal grid points per period. This should be an odd number.')
    parser.add_argument('--NzetaScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nzeta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Ntheta', type=int, nargs=1, required=False, default=[43], help='Number of poloidal grid points. This should be an odd number.')
    parser.add_argument('--NthetaScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Ntheta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nxi', type=int, nargs=1, required=False, default=[125], help='Number of Legendre polynomials used to represent the pitch-angle dependence of the distribution function.')
    parser.add_argument('--NxiScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nxi that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nx', type=int, nargs=1, required=False, default=[7], help='Number of grid points in energy used to represent the distribution function.')
    parser.add_argument('--NxScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nx that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--NL', type=int, nargs=1, required=False, default=[4], help='Number of Legendre polynomials used to represent the Rosenbluth potentials. Increasing this hardly changes the results, so it can almost certainly be left alone.')
    parser.add_argument('--NLScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of NL that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--solverTol', type=float, nargs=1, required=False, default=[1e-6], help='Tolerance used to define convergence of the iterative (Krylov) solver.')
    parser.add_argument('--solverTolScan', type=float, nargs=2, required=False, default=[0.1, 10.0], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of solverTolerance that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--saveLoc', type=str, nargs='*', required=False, default=[None], help='Location(s) in which to save written files - this will act as the main directory(ies) for a set of SFINCS runs. Defaults to either <profilesIn> or <eqIn> location(s) -- whichever input has more locations will be chosen as the default. If len(<profilesIn>) == len(<eqIn>), defaults to <profilesIn> location(s). If you input multiple files, order matters! Note that if you specify 1 <saveLoc> and multiple <profilesIn> or <eqIn>, the code will attempt to save all the generated files in the same directory. Due to the current (strict) naming conventions of SFINCS, this is probably not useful because the last-written files will overwrite their predecessors, but the feature is included for completeness.')
    parser.add_argument('--nNodes', type=int, nargs=1, required=False, default=[None], help='Total number of nodes to use for each SFINCS run. If you do not use <noRun>, you must specify at least one of <nNodes> and <nTasks>.')
    parser.add_argument('--nTasksPerNode', type=int, nargs=1, required=False, default=[None], help='Number of MPI tasks to use on each node for each SFINCS run. This parameter should only be used if <nNodes> is specified and should not be used with <nTasks>.')
    parser.add_argument('--nTasks', type=int, nargs=1, required=False, default=[None], help='Total number of MPI tasks to use for each SFINCS run. If you do not use <noRun>, you must specify at least one of <nNodes> and <nTasks>.')
    parser.add_argument('--mem', type=int, nargs=1, required=False, default=[None], help='Total amount of memory (MB) allocated for each SFINCS run.')
    parser.add_argument('--time', type=str, nargs=1, required=False, default=['00-06:00:00'], help='Wall clock time limit for the batch runs. Format is DD-HH:MM:SS. Note that SFINCS typically has the most trouble converging near the magnetic axis (due to the lower collisionality there cause by peaked temperature profiles), so you may need to increase <time> for runs near the axis.')
    parser.add_argument('--noProfiles', action='store_true', default=False, help='Do not write a profiles file.')
    parser.add_argument('--noNamelist', action='store_true', default=False, help='Do not write an input.namelist file.')
    parser.add_argument('--noBatch', action='store_true', default=False, help='Do not write a job.sfincsScan file.')
    parser.add_argument('--noRun', action='store_true', default=False, help='Do not run sfincsScan.')
    parser.add_argument('--notifs', type=str, nargs=1, required=False, default=['bad'], help='Dictate which Slurm notification emails you would like to receive. By default, you will only receive emails when something bad happens to your job (such as a failure). You may also specify "all" or "none", which have the (intuitive) meanings indicated in the Slurm documentation. Note that the environment variable SFINCS_BATCH_EMAIL must be set for <notifs> to work correctly.')
    parser.add_argument('--noConfirm', action='store_true', default=False, help='Instruct sfincsScan to create folders and jobs without asking for confirmation first.')
    args = parser.parse_args()

    if not all([i in ['sym', 'asym'] for i in args.bcSymmetry]):
        raise IOError('Each element of <bcSymmetry> must be set to either "sym" or "asym".')

    if len(args.bcSymmetry) not in [1, len(args.eqIn)]:
        raise IOError('The length of <bcSymmetry> must either be 1 or the same as <eqIn>.')

    if args.loadPot:
        if args.numErSubscan[0] < 1:
            raise IOError('If you activate <loadPot>, you must have <numErSubscan> >= 1.')
        if args.radialGradientVar[0] != 1:
            raise IOError('If you activate <loadPot>, you must specify <radialGradientVar> = 1.')

    if args.minSeedEr[0] > args.maxSeedEr[0]:
        raise IOError('<minSeedEr> must be less than or equal to <maxSeedEr>.')

    if args.minSeedEr[0] == 0 and args.maxSeedEr[0] == 0 and args.numErSubscan[0] > 1:
        raise IOError('If you set <minSeedEr> = <maxSeedEr> = 0 and numErSubscan > 1, you will create redundant runs.')

    if not args.loadPot and (args.minSolverEr[0] > args.minSeedEr[0] or args.maxSolverEr[0] < args.maxSeedEr[0]):
        errStr = 'The range of Er available to ambipolarSolve (set by <*SolverEr>) is likely smaller than the range '
        errStr += 'you wish to explore (set by <*SeedEr>). '
        errStr += 'The comparison made to throw this error is approximate because the units of <*SolverEr> are fixed while '
        errStr += 'those of <*SeedEr> are variable, but there should still be plenty of "space" around the <*SeedEr> range.'
        raise IOError(errStr)

    if not args.loadPot and not (args.minSolverEr[0] < args.seedEr[0] < args.maxSolverEr[0]):
        errStr = 'It appears that <seedEr> was set outside the range specified by the <*SolverEr> pair. '
        errStr += 'The comparison made to throw this error is approximate because the units of <*SolverEr> are fixed while '
        errStr += 'those of <seedEr> are variable, but there should still be plenty of "space" around <seedEr>.'
        raise IOError(errStr)
    
    if args.Nyquist[0] not in [1,2]:
        raise IOError('An invalid <Nyquist> choice was specified. Valid inputs are the integers 1 and 2.')

    if args.radialVar[0] not in [0,1,2,3]:
        raise IOError('An invalid <radialVar> choice was specified. Valid inputs are the integers 0, 1, 2, and 3.')
    
    if args.radialGradientVar[0] not in [0,1,2,3,4]:
        raise IOError('An invalid <radialGradientVar> choice was specified. Valid inputs are the integers 0, 1, 2, 3, and 4.')

    if args.Nzeta[0]%2 == 0:
        raise IOError('<Nzeta> should be odd.')
    
    if args.Ntheta[0]%2 == 0:
        raise IOError('<Ntheta> should be odd.')
    
    lens = [len(args.profilesIn), len(args.eqIn), len(args.saveLoc)]
    maxLen = max(lens)
    for length in lens:
        if length != 1 and length != maxLen:
            raise IOError('Regarding <profilesIn>, <eqIn>, and <saveLoc>: any of these three inputs with length greater than 1 must have the same length as the other inputs with length greater than 1.')

    if not args.noRun:

        if args.nNodes[0] is None and args.nTasks[0] is None:
            raise IOError('You must specify at least one of <nNodes> and <nTasks>.')

        if args.nNodes[0] is None and args.nTasksPerNode[0] is not None:
            raise IOError('You cannot specify <nTasksPerNode> without also specifying <nNodes>.')

        if args.nTasks[0] is not None and args.nTasksPerNode[0] is not None:
            raise IOError('You cannot specify both <nTasks> and <nTasksPerNode>.')

    strippedTime = args.time[0].strip()
    timeDaySplit = strippedTime.split('-')
    
    assert len(timeDaySplit) == 2, 'There appears to be more than one "-" character in the <time> input.'

    try:
        int(timeDaySplit[0])
    except ValueError:
        raise IOError('It appears that the "days" portion of the <time> input is incorrect.')

    timeHourSplit = timeDaySplit[-1].split(':')

    try:
        _ = [int(elem) for elem in timeHourSplit]
    except ValueError:
        raise IOError('It appears that at least one of the "hours", "minutes", or "seconds" portions of the <time> input is incorrect.')

    if args.notifs[0].lower() not in ['bad', 'all', 'none']:
        raise IOError('Invalid option specified for <notifs>. Valid options are "bad", "all", and "none".')

    return args

def getPlotArgs():

    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for plotting SFINCS outputs.
    '''

    import argparse
    from os.path import isdir
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sfincsDir', type=str, nargs='*', required=True, help='Top directory(ies) for SFINCS run(s), with path(s) if necessary. Such directories contain subdirectories which either contain SFINCS output files (*.h5) or more subdirectories for the electric field scan. In the latter case, those subsubdirectories contain SFINCS output files. Note that the "most complete" distribution function available ("vm" for calculations without Phi1 and "vd" for calculations with Phi1) will be used for most plots. If you input multiple directories, order matters!')
    parser.add_argument('--radialVar', type=int, nargs=1, required=False, default=[3], help='ID of the radial coordinate used in the input.namelist file to specify which surfaces should be scanned over. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, and 3 = rN (which is the STELLOPT rho)')
    parser.add_argument('--radialVarBounds', type=float, nargs=2, required=False, default=[-1, -1], help='Two floats, which are (in order) the minimum and maximum values of <radialVar> that will be plotted. If one of the inputs is negative, it will be ignored (so that the min or max is not limited).')
    parser.add_argument('--saveLoc', type=str, nargs='*', required=False, default=[None], help='Location(s) in which to save plots, plot data, and informational *.txt files. Defaults to <sfincsDir>/processed/. If you input multiple directories, order matters!')
    parser.add_argument('--checkConv', action='store_true', default=False, help='Instead of plotting anything, just check if the SFINCS runs in the <sfincsDir> location(s) converged. If they all did, you will receive no output.')
    args = parser.parse_args()

    if not all([isdir(item) for item in args.sfincsDir]):
        raise IOError('The inputs given in <sfincsDir> must be directories.')

    if args.radialVar[0] not in [0,1,2,3]:
        raise IOError('An invalid <radialVar> choice was specified. Valid inputs are the integers 0, 1, 2, and 3.')
    
    lens = [len(args.sfincsDir), len(args.saveLoc)]
    maxLen = max(lens)
    for length in lens:
        if length != 1 and length != maxLen:
            raise IOError('If both <sfincsDir> and <saveLoc> have length greater than 1, they must be the same length.')
    
    return args

def getPhi1SetupArgs():
    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for setting up SFINCS runs with Phi1 included.
    '''
    
    import argparse
    from os.path import isdir

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sfincsDir', type=str, nargs='*', required=True, help='Top directory(ies) for SFINCS run(s), with path(s) if necessary. Each directory must contain radial (or radial and electric field) subdirectories, each with an "input.namelist" file, "job.sfincsScan" file, and "sfincsOutput.h5" file. These files will simply be copied and modified as necessary to include Phi1 calculations. If you input multiple directories, order matters!')
    parser.add_argument('--saveLoc', type=str, nargs='*', required=False, default=[None], help='Top-level directory(ies) in which to save modified files and informational *.txt files. The directory structure will be copied from <sfincsDir>. Defaults to <sfincsDir>+"_Phi1". If you input multiple directories, order matters!')
    parser.add_argument('--noRun', action='store_true', default=False, help='Copy/write files, but do not launch SFINCS.')
    args = parser.parse_args()

    if not all([isdir(item) for item in args.sfincsDir]):
        raise IOError('The inputs given in <sfincsDir> must be directories.')
    
    lens = [len(args.sfincsDir), len(args.saveLoc)]
    maxLen = max(lens)
    for length in lens:
        if length != 1 and length != maxLen:
            raise IOError('If both <sfincsDir> and <saveLoc> have length greater than 1, they must be the same length.')
    
    return args

def getAxisParamsArgs():
    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for getting axis parameters from VMEC wout files.
    '''

    import argparse
    from os.path import isdir

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wout', type=str, nargs=1, required=True, help='wout file from which to pull axis information. Note that stellarator symmetry is assumed!')
    args = parser.parse_args()

    if isdir(args.wout[0]):
        raise IOError('The input given in <wout> must be a file.')

    return args

def getCompoundPlotArgs():
    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for creating compound plots.
    '''
    
    import argparse
    from os.path import isdir

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--data', nargs='+', required=True, help='*.dat files to be plotted, with addresses if necessary. The first column of each file must be horizontal coordinate values, while the remaining columns must be corresponding vertical coordinate (data) values. The file plot.py produces properly-structured data files automatically. Note that you should choose *.dat files with the same horizontal coordinate. Each column (except the first) in the first file passed to this argument will become a curve in the output plot, then each column (except the first) in the second file, and so on -- this is how one can specify the order of the <legend> and <colors> arguments, for example.')
    parser.add_argument('--plotType', type=str, nargs=1, required=False, default=['linear'], help='Type of vertical axis. Options are "linear" and "semilogy".')
    parser.add_argument('--xlabel', type=str, nargs=1, required=False, default=[''], help='Label for horizontal axis. Be sure to write in quotes!')
    parser.add_argument('--ylabel', type=str, nargs=1, required=False, default=[''], help='Label for vertical axis. Be sure to write in quotes!')
    parser.add_argument('--xtick', type=float, nargs=1, required=False, default=[0.1], help='Sets spacing of the horizontal ticks.')
    parser.add_argument('--xmin', type=float, nargs=1, required=False, default=[None], help='Sets minimum value on the horizontal axis.')
    parser.add_argument('--legend', nargs='+', required=False, default=[''], help='Legend entries. Can accept one or more arguments. Be sure to write each argument in quotes!')
    parser.add_argument('--saveLoc', type=str, required=False, default='.', help='Location in which to save the produced plot. Defaults to the current working directory.')
    parser.add_argument('--fileName', type=str, required=False, default='compositePlot', help='Name of the produced plot without a file extension.')
    parser.add_argument('--fileType', type=str, required=False, default='pdf', help='Type of file you wish to save. Options are "pdf" and "png".')
    parser.add_argument('--fontSize', type=float, required=False, default=18, help='Font size for plots.')
    parser.add_argument('--ymargin', type=float, required=False, default=None, help='Set vertical axis autoscaling margin. Must be between 0 and 1.')
    parser.add_argument('--colors', nargs='+', required=False, default=[''], help='Curve colors. For possible options, see the matplotlib documentation. Note that the default matplotlib colors are from the Tableau palette -- for instance, the default blue is "tab:blue", the default orange is "tab:orange", and so forth. Can accept one or more arguments. Be sure to write each argument in quotes! If one argument is input, all curves will have the same color. If multiple arguments are input, this must be the same as the number of lines that are to be plotted.')
    parser.add_argument('--lineStyles', nargs='+', required=False, default=[''], help='Curve line styles. For possible options, see the matplotlib documentation. It is suggested that you use the spelled-out forms for the line styles rather than their abbreviations (e.g. use "solid", "dashed", etc.). Can accept one or more arguments. Be sure to write each argument in quotes! If one argument is input, all curves will have the same line style. If multiple arguments are input, this must be the same as the number of lines that are to be plotted.')
    parser.add_argument('--markers', nargs='+', required=False, default=[''], help='Curve marker styles. For possible options, see the matplotlib documentation. Can accept one or more arguments. Be sure to write each argument in quotes! If one argument is input, all curves will have the same marker style. If multiple arguments are input, this must be the same as the number of lines that are to be plotted.')
    args = parser.parse_args()
    
    for item in args.data:
        if isdir(item):
            raise IOError('The input(s) given in <data> must be files.')

    if args.plotType[0] not in ['linear', 'semilogy']:
        raise IOError('<plotType> must be either "linear" or "semilogy"')

    if args.fileType not in ['pdf', 'png']:
        raise IOError('<fileType> must be either "pdf" or "png"')

    if args.ymargin is not None and not (0 <= args.ymargin <= 1):
        raise IOError('<ymargin> must be between 0 and 1.')

    return args
    
def getChooseErArgs():

    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts for choosing the right radial electric field on a given flux surface.
    '''

    import argparse
    from os.path import isdir
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sfincsDir', type=str, nargs=1, required=True, help='Top directory for SFINCS run, with path if necessary. This directory must contain subdirectories which contain more subdirectories for the electric field scan. The lowest-level subsubdirectories contain SFINCS output files (*.h5).')
    parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=[None], help='Location in which to save outputs. If <filter> is not used, these outputs are plots and informational *.txt files. In this case, the default is <sfincsDir>/determineEr/. If <filter> is used, the output is a "mirror" of <sfincsDir> containing only the "correct" electric field information. In this case, the default is <sfincsDir>+"_correctEr". The default is recommended, particularly when <filter> is not used.')
    parser.add_argument('--filter', action='store_true', default=False, help='Once all the roots for all flux surfaces of interest in a given <sfincsDir> are determined, this option can be used to copy only the subdirectories that contain the "correct" electric field information from <sfincsDir> to <saveLoc>. If plot.py is then run on <saveLoc>, the "true" behavior of the system will be seen. Note that <sfincsDir> must contain a determineEr/ subdirectory with a rootsToUse.txt file for this option to work.')
    parser.add_argument('--noAmbiSolve', action='store_true', default=False, help='Disable ambipolarSolve. This means that guessed value(s) of Er will be used without modification during the SFINCS run(s). Note that if Phi1 effects are included in the calculations, you MUST toggle <noAmbiSolve> and repeatedly use this script to find the roots.')
    parser.add_argument('--noRun', action='store_true', default=False, help='Copy files, but do not launch SFINCS runs.')
    parser.add_argument('--maxRootJr', type=float, nargs=1, required=False, default=[1.0e-12], help='Maximum radial current (defined as in SFINCS) that may be present for a given electric field value to be considered a "root". The default is the same as that used in writeNamelist.py and is recommended.')
    args = parser.parse_args()

    if not isdir(args.sfincsDir[0]):
        raise IOError('The input given in <sfincsDir> must be a directory.')
    
    return args

def getFileInfo(inFile, saveLoc, outFileName):

    '''
    Inputs:
        inFile: String with (relative or absolute) path
                to an input file.
        saveLoc: String with (relative or absolute) path 
                 where other files (such as the outputs
                 of other scripts) should be saved. Defaults
                 to the location of inFile.
        outFileName: String with name for an outFile.
    Outputs:
        Strings with the inFile absolute path, inFile name,
        inFile path, outFile path, and outFile absolute path.
    '''

    from os.path import abspath, dirname, basename, join

    inFile = abspath(inFile)
    inFilePath = dirname(inFile)
    inFileName = basename(inFile)

    if saveLoc == None:
        outFilePath = inFilePath
    else:
        outFilePath = abspath(saveLoc)

    outFile = join(outFilePath, outFileName)

    return inFile, inFileName, inFilePath, outFilePath, outFile

def cleanStrings(inputList):

    '''
    Inputs:
        List of strings.
    Outputs:
        inputList, but without spaces and with all other 
        characters lowercase.
    '''

    outList = []
    for item in inputList:
        cleaned = item.strip().lower()
        outList.append(cleaned)

    return outList

def listifyBEAMS3DFile(inputFile):
    
    '''
    Inputs:  
        inputFile: STELLOPT input.namelist file with a BEAMS3D section.
    Outputs: 
        The variable assignment data in the BEAMS3D section
        as a list of lists. Each sublist contains the variable
        name as its first element. Each subsequent element is 
        a value assigned to the variable. This means that 
        singly-valued variables will have lists with length 2
        and vector-valued variables will have lists with length 
        len(vector)+1. Note that any duplicate sublists are
        filtered such that only one copy remains.
    '''
    
    from itertools import groupby
        
    with open(inputFile,'r') as f:
        beams3dSectionStartFlag = False
        dataLines = []
        
        for line in f:
            cleaned = line.strip().lower()
            if (beams3dSectionStartFlag == True):
                if cleaned == '/':
                    break
                dataLines.append(cleaned)
            if cleaned == '&beams3d_input':
                beams3dSectionStartFlag = True
    
    listifiedData = []

    for line in dataLines:
        
        if line[0] == '!':
            continue
        
        precleaned = [i.strip() for i in line.split('=')]

        if len(precleaned) != 2:
            raise IOError('The script thinks that there were two equals signs in a variable assignment line in {}. Something is wrong.'.format(inputFile))
        
        if (' ' or '\t') in precleaned[1]:
            cleaned = [precleaned[0]] + precleaned[1].split()
        else:
            cleaned = precleaned

        listifiedData.append(cleaned)

    listifiedData.sort()

    redundanciesRemoved = list(listifiedData for listifiedData,_ in groupby(listifiedData))

    return redundanciesRemoved

def makeProfileNames(listOfPrefixes):

    '''
    Inputs:
        listOfPrefixes: A list of the form ['name1','name2',...].
                        Typical names are NE, TI, and so forth.
                        Note that repeated inputs are possible -
                        this is useful for, say, using the same
                        profiles for two different variables,
                        such as NI and NE.
    Outputs:
        A list of lists containing BEAMS3D variables derived 
        from listOfPrefixes, such as
        [[name1_AUX_S, name1_AUX_F], ...].
    '''
  
    output_names = []
    for prefix in listOfPrefixes:
        s = prefix + '_aux_s'
        f = prefix + '_aux_f'
        appendor = [s, f]
        output_names.append(appendor)

    return output_names

def extractProfileData(dataList, nameList):

    '''
    Inputs:  
        dataList: A list of lists, as from the listifyBEAMS3DFile function,
                   with a string as the first element and length >= 2.
        nameList: A list of lists, as from the makeProfileNames function. 
                  Each sublist contains a pair of strings to search for in 
                  dataList.
    Outputs:
        A dictionary. Each key is a unique prefix from nameList. Each value
        contains a dictionary with two entries. The 'iv' entry is a list
        of lists with the radial coordinate (normalized toroidal flux) for
        the given variable. The 'dv' entry is a lists of lists with values
        corresponding to those radial coordinates. If only a single species
        (or Phi) is specified for a given prefix, the lists will look like
        [[...values...]]. With multiple species, they will look like
        [[...values1...],[...values2...],...]. Note that the dimensions of
        the 'iv' and 'dv' arrays will always match. This means that some
        'iv' values may be repeated, such as when multiple ion profiles
        are specified using the same set of radial coordinates.
    '''

    import warnings

    matched = []
    dataDict = {}
    for namePair in nameList:
        
        strippedName = namePair[0].split('_')[0]

        matchedPair = {}
        for name in namePair:
            foundMatch = False

            allSpeciesData = {}
            for dataVec in dataList:

                # We need to isolate the name of the variable
                dataLabel = dataVec[0]
                splitAtOpenPar = dataLabel.split('(')
                dataName = splitAtOpenPar[0]
                
                if dataName == name:
                    # We need to sort out the variable indices
                    if len(splitAtOpenPar) == 1: # The label has no indices, so its array has only one row
                        speciesIndex = 0
                    
                    elif len(splitAtOpenPar) == 2: # The label has indices that we need to sort out
                        indices = [item.strip() for item in splitAtOpenPar[1].split(')')[0].split(',')]
                        speciesIndex = int(indices[0]) - 1 # Python indices start from 0, STELLOPT indices start from 1
                        if indices[1] != ':':
                            raise IOError('Piecewise indexing for profile declarations (as in {}) is not yet supported.'.format(dataLabel))
                    
                    else:
                        raise IOError('In the variable declaration line for {}, the "(" character apparently appears twice. Something is wrong.'.format(dataLabel))

                    # Now assign the data as an IV or a DV
                    allSpeciesData[speciesIndex] = [float(i) for i in dataVec[1:]]
                    foundMatch = True
            
            allSpeciesDataList = list(dict(sorted(allSpeciesData.items())).values()) # Note that we will always have a list of lists
            if name[-1] == 's':
                matchedPair['iv'] = allSpeciesDataList
            elif name[-1] == 'f':
                matchedPair['dv'] = allSpeciesDataList
            else:
                raise IOError('The read variable suffix for {} is not "S" or "F". Something is wrong.'.format(dataLabel))

            if not foundMatch:
                warnings.warn('No match could be found for the variable "{}" in the given dataList!'.format(name))

        # Clean up the empty lists in matchedPair
        cleanMatchedPair = {}
        for key,data in matchedPair.items():
            cleanMatchedPair[key] = [item for item in data if item != []]
 
        # Store the data
        matched.append(cleanMatchedPair)
        dataDict[strippedName] = cleanMatchedPair

    # Make the IV and DV parts of cleanMatchedPair the same length
    for key,data in dataDict.items():
        ivLen = len(data['iv'])
        dvLen = len(data['dv'])
        if ivLen == dvLen:
            pass
        elif ivLen == 1 and dvLen > 1:
            dataDict[key]['iv'] = data['iv'] * dvLen
        elif ivLen > 1 and dvLen == 1:
            raise IOError('Multiple independent variable sets were specified for one dependent variable set in the {} array. Something is wrong.'.format(strippedName))
        else:
            raise IOError('The structure of the data in the {} array is irregular. Something is wrong.'.format(strippedName))
    
    if not any(matched):
        raise IOError('No searched variables were found.')
    
    return dataDict

def extractScalarData(dataList, nameList):

    '''
    Inputs:  
        dataList: A list of lists, as from the listifyBEAMS3DFile function,
                   with a string as the first element and length >= 2.
        nameList: A list of variable names to search for in dataList. 
    Outputs:
        A dictionary. Each key is a brief name of a variable from nameList.
        Each value contains a list with one or more floats corresponding to
        the variable.
    '''

    import warnings
    
    matched = []
    dataDict = {}
    for name in nameList:
        
        foundMatch = False
        strippedName = name.split('_')[-1]
        
        for dataVec in dataList:
            
            if dataVec[0] == name:
                floats = [float(i) for i in dataVec[1:]]
                dataDict[strippedName] = floats
                foundMatch = True
                matched.append(floats)
                break

        if not foundMatch:
            warnings.warn('No match could be found for the variable "{}" in the given dataList!'.format(name))
        
    if not any(matched):
        raise IOError('No searched variables were found.')
    
    return dataDict

def sortProfileFunctions(inDict):

    '''
    Inputs:
        A dictionary (such as from the interpolatedData function) whose keys include 'ne', 'ni', 'te', and 'ti',
        and whose values are lists of functions corresponding to the profiles of each species. It is expected
        that the values of the *e variables have length 1 since electrons are a unique species. The lengths
        of the 'ni' and 'ti' values can be arbitrary (since one can include as many ion species in the calculation as
        they wish), subject to the constraint that both are either the same length (unique profiles for each species)
        or one is length 1 (the same profile will be used for all species).
    Outputs:
        A list with the interpolation functions ordered in the form [N1, T1, N2, T2, ...] suitable for use in
        the generateDataText function. Note that electrons are always species 1! 
    '''

    if len(inDict['ne']) != 1 or len(inDict['te']) != 1 :
        raise IOError('The lengths of the "ne" and "te" values must be 1.')

    lenNI = len(inDict['ni'])
    lenTI = len(inDict['ti'])

    if (lenNI != lenTI) and (lenNI != 1 and lenTI != 1):
        raise IOError('The lengths of "ni" and "ti" must either be the same, or one of them must have length 1.')

    outList = [inDict['ne'][0], inDict['te'][0]]

    if lenNI > lenTI:
        NIs = inDict['ni']
        TIs = inDict['ti'] * lenNI
    elif lenNI < lenTI:
        NIs = inDict['ni'] * lenTI
        TIs = inDict['ti']
    else:
        NIs = inDict['ni']
        TIs = inDict['ti']

    for (NI, TI) in zip(NIs, TIs):
        outList.append(NI)
        outList.append(TI)

    return outList

def generatePreamble(radial_coordinate_ID=1):

    '''
    Inputs:
        The (integer) ID for the radial coordinate used in profiles file:
        0 = psiHat, 1 = psiN, 2 = rHat, 3 = rN.
        Note that BEAMS3D always uses the normalized toroidal flux, which
        corresponds to radial_coordinate_ID=1.
    Outputs:
        String containing the preamble to the profiles file.
    '''

    stringToWrite = '# This is an integer specifying the radial coordinate used in this file, which can be different from the one specified by inputRadialCoordinate in input.namelist (<radialVar>).\n'
    stringToWrite += '{}\n'.format(radial_coordinate_ID)
    stringToWrite += '# The following lines contain profile information in this format:\n'
    stringToWrite += '# radius\tNErs\tgeneralEr_min\tgeneralEr_max\tnHat(species 1)\tTHat(species 1)\tnHat(species 2)\tTHat(species 2)\t...\n'
    stringToWrite += '# The meaning of the generalEr_* variables is set by inputRadialCoordinateForGradients in input.namelist (<radialGradientVar>).\n'

    return stringToWrite

def generateDataText(radii, *funcs):

    '''
    Inputs:
        radii: A list of radii at which to evaulate *funcs.
        *funcs: Functions that can take one element of radii
                at a time and output a single value (each) 
                needed in the profiles file.
    Outputs:
        Text constituting all the computer-readable information 
        in profiles file except for radial_coordinate_ID.
    '''
    
    def stringifyItem(item, endOfLine=False):
        
        if endOfLine:
            out = str(item) + '\n' # Note that the last line in the file should have a \n (UNIX standard)
        else:
            out = str(item) + '\t'

        return out
        
    stringOut = ''
    for radius in radii:
        stringOut += stringifyItem(radius)
        for func in funcs[:-1]:
            stringOut += stringifyItem(func(radius))
        stringOut += stringifyItem(funcs[-1](radius), endOfLine=True)

    return stringOut

def writeFile(outFile, stringToWrite, silent=False):

    '''
    Inputs:
        outFile: absolute name (with path) of file
                 to be written.
        stringToWrite: string to write in outFile.
        silent: if True, suppresses the output message
                that outFile has been written.
    Outputs:
        [A written file, and possibly a notification message.]
    '''
   
    _, outFileName, _, _, _ = getFileInfo(outFile, '/arbitrary/path/', 'arbitrary')

    with open(outFile, 'w') as f:
        f.write(stringToWrite)

    if not silent:
        messagePrinter('{} file written.'.format(outFileName))

def findFiles(name, path, raiseError=False):

    '''
    Inputs:
        name: string with the name of a file to search for.
              Note that the name must be exact (there is
              no name matching).
        path: string with path of directory to search 
              (recursively) for name.
        raiseError: if True, raise an IOError when no files
                    called name are found in path
    Outputs:
        Sorted list with absolute paths to files called
        name within path.
    '''
    
    from os import walk
    from os.path import join

    result = []
    for root, dirs, files in walk(path):
        if name in files:
            result.append(join(root, name))
    
    result.sort() # Not necessary, just makes outputs a bit easier to follow
    
    if raiseError and len(result) == 0:
        raise IOError('No files called {} could be found in {}.'.format(name, path))
    
    return result

def adjustInputLengths(inListDict):

    '''
    Inputs:
        A dictionary containing lists. Any of these lists
        with length greater than 1 must have the same length
        as the other lists with length greater than 1.
    Outputs:
        outListDict: the same dictionary as inListDict, but
                     with each list being the same length.
        longLists: a list with the keys of inListDict that
                   had the largest length.
        maxLen: length of longest list in inListDict, and
                therefore the length of all lists in
                outListDict.
    '''

    maxLen = max([len(data) for key,data in inListDict.items()])
    
    outListDict = {}
    longLists = []
    for key,data in inListDict.items():
        if len(data) != maxLen:
            outListDict[key] = data * maxLen
        else:
            outListDict[key] = data
            longLists.append(key)

    return outListDict, longLists, maxLen

def makeDir(saveLoc):

    '''
    Inputs:
        Relative or absolute path to desired save location.
    Outputs:
        [Desired save location is created if not already present.]
        The absolute address of the created directory is output.
    '''

    from os import makedirs

    _, _, _, outDir, _ = getFileInfo('/arbitrary/path', saveLoc, 'arbitrary')
    makedirs(outDir, exist_ok=True)

    return outDir

def radialVarDict():

    '''
    Inputs:
        [None]
    Outputs:
        Dictionary with the indices of the SFINCS radial variables
        as keys and the written-out forms of the variables as values.
    '''

    return {0:'psiHat', 1:'psiN', 2:'rHat', 3:'rN', 4:'rHat'}

def prettyRadialVar(inString, innerOnly=False):
    
    '''
    Inputs:
        inString: string containing a value from radialVarDict.
        innerOnly: if False, all the 'extras' needed for the
                   output to act as, say, an axis label will
                   be included. If True, these 'extras' will
                   not be included (which is useful if the
                   output is to be included in, say, another
                   equation).
    Outputs:
        Raw string that can be used to print a formatted
        form of inString, such as with Matplotlib.
    '''

    if inString == 'psiHat':
        core = r'\hat{\psi}'
        if innerOnly:
            return core
        else:
            return r'${}$'.format(core)
    elif inString == 'psiN':
        core = r'\psi_{N}'
        if innerOnly:
            return core
        else:
            return r'$%s = s = \left( r_{\mathrm{eff}}/a_{\mathrm{eff}}\right)^{2}$'%core
    elif inString == 'rHat':
        core = r'\hat{r}'
        if innerOnly:
            return core
        else:
            return r'${}$'.format(core)
    elif inString == 'rN':
        core = r'r_{N}'
        if innerOnly:
            return core
        else:
            return r'$%s = \rho = r_{\mathrm{eff}}/a_{\mathrm{eff}}$'%core
    else:
        raise IOError('Invalid radial coordinate input.')

def prettyDataLabel(inString):

    '''
    Inputs:
        String containing a variable name from the SFINCS
        output (*.h5) file. Note that only a few variables
        are currently supported.
    Outputs:
        Raw string that can be used to print a formatted
        form of inString, such as with Matplotlib.
    '''

    if '_' not in inString:
        
        PhiHat = r'$\hat{\Phi}$'
        def derFormat(top, bottom):
            return r'$\frac{d %s}{d %s}$' % (top, bottom)
        
        extensiveParticleFluxUnits = r' $\mathrm{\left(\frac{1}{s}\right)}$'
        extensiveHeatFluxUnits = r' $\mathrm{\left(W\right)}$'
        extensiveMomentumFluxUnits = r' $\mathrm{\left(\frac{kg T m}{s^{-2}}\right)}$'
        extensiveRadialCurrentUnits = r' $\mathrm{\left(A\right)}$'

        if inString == 'Er':
            return r'Radial electric field $\mathrm{\left(\frac{V}{m}\right)}$'

        elif inString == 'dPhiHatdpsiHat':
            return 'Radial electric field ' + derFormat(PhiHat, prettyRadialVar('psiHat', innerOnly=True))

        elif inString == 'dPhiHatdpsiN':
            return 'Radial electric field ' + derFormat(PhiHat, prettyRadialVar('psiN', innerOnly=True))
        
        elif inString == 'dPhiHatdrHat':
            return 'Radial electric field ' + derFormat(PhiHat, prettyRadialVar('rHat', innerOnly=True))
        
        elif inString == 'dPhiHatdrN':
            return 'Radial electric field ' + derFormat(PhiHat, prettyRadialVar('rN', innerOnly=True))
        
        elif inString == 'FSABFlow':
            return r'FSAB parallel flow $\mathrm{\left(\frac{T}{m^{2} s}\right)}$'
        
        elif inString == 'FSABjHat':
            return r'FSAB bootstrap current $\mathrm{\left(\frac{T A}{m^{2}}\right)}$'

        elif inString in ['FSABjHatOverRootFSAB2', 'FSABjHatOverB0']:
            return r'Bootstrap current $\mathrm{\left(\frac{A}{m^{2}}\right)}$'
        
        elif inString == 'extensiveParticleFlux':
            return r'Neoclassical particle flux' + extensiveParticleFluxUnits

        elif inString == 'extensiveClassicalParticleFlux':
            return r'Classical particle flux' + extensiveParticleFluxUnits
        
        elif inString == 'extensiveTotalParticleFlux':
            return r'Total particle flux' + extensiveParticleFluxUnits
        
        elif inString == 'extensiveHeatFlux':
            return r'Neoclassical energy flux' + extensiveHeatFluxUnits
        
        elif inString == 'extensiveClassicalHeatFlux':
            return r'Classical energy flux' + extensiveHeatFluxUnits
        
        elif inString == 'extensiveTotalHeatFlux':
            return r'Total energy flux' + extensiveHeatFluxUnits
        
        elif inString == 'extensiveMomentumFlux':
            return r'Neoclassical momentum flux' + extensiveMomentumFluxUnits
        
        elif inString == 'extensiveRadialCurrent':
            return r'Radial current' + extensiveRadialCurrentUnits
        
        else:
            raise IOError('Formatting has not yet been specified for the variable {}.'.format(inString))
    
    else:
        
        # Sort out the parts you need to make sense of inString and do some basic administrative checks
        parts = inString.split('_')
        
        if len(parts) not in [2, 3]:
            raise IOError('Formatting has not yet been specified for the variable {}.'.format(inString))

        if len(parts) == 3 and parts[1] != 'vm' and parts[1] != 'vd':
            raise IOError('Formatting has not yet been specified for the variable {}.'.format(inString))
        
        label = parts[0]
        radVar = prettyRadialVar(parts[-1], innerOnly=True)
        
        # Write some label strings that will be used below
        directionStatement = r' in $\nabla {}$ direction '.format(radVar) # Note that this includes spaces on either side for convenience
        particleFluxUnits = r'$\mathrm{\left(\frac{1}{m^{3} s}\right)}$'
        heatFluxUnits = r'$\mathrm{\left(\frac{W}{m^{3}}\right)}$'
        momentumFluxUnits = r'$\mathrm{\left(\frac{kg T}{m^{2} s^{2}}\right)}$'
        radialCurrentUnits = r'$\mathrm{\left(\frac{A}{m^{3}}\right)}$'

        # Write the output
        if label == 'particleFlux':
            return r'Neoclassical particle flux' + directionStatement + particleFluxUnits

        elif label == 'classicalParticleFlux':
            return r'Classical particle flux' + directionStatement + particleFluxUnits
        
        elif label == 'classicalParticleFluxNoPhi1':
            return r'Classical particle flux (neglecting $\Phi_{1}$)' + directionStatement + particleFluxUnits

        elif label == 'totalParticleFlux':
            return r'Total particle flux' + directionStatement + particleFluxUnits

        elif label == 'heatFlux':
            return r'Neoclassical energy flux' + directionStatement + heatFluxUnits

        elif label == 'classicalHeatFlux':
            return r'Classical energy flux' + directionStatement + heatFluxUnits
        
        elif label == 'classicalHeatFluxNoPhi1':
            return r'Classical energy flux (neglecting $\Phi_{1}$)' + directionStatement + heatFluxUnits

        elif label == 'totalHeatFlux':
            return r'Total energy flux' + directionStatement + heatFluxUnits

        elif label == 'momentumFlux':
            return r'Neoclassical momentum flux' + directionStatement + momentumFluxUnits

        elif label == 'radialCurrent':
            return r'Radial current' + directionStatement + radialCurrentUnits

        else:
            raise IOError('Formatting has not yet been specified for the variable {}.'.format(inString))

def messagePrinter(inString):

    '''
    Inputs:
        A string to be formatted in a consistent fashion
        and printed. This is primarily intended for use
        with informational statements.
    Outputs:
        [Printing a formatted version of inString.]
    '''

    print('***StelloptPlusSfincs: %s***'%inString)

def now():
    
    '''
    Inputs:
        [None]
    Outputs:
        A string containing the current time.
    '''

    import datetime

    t = lambda: str(datetime.datetime.now())
    
    return t()
    
def saveTimeStampFile(saveLoc, fileName, initialString):
    
    '''
    Inputs:
        saveLoc: The absolute directory in which a file should
                 be written.
        fileName: The name of the file to be written, excluding
                  the ".txt" extension.
        initialString: A string that will be included before
                       the actual time stamp.

    Outputs:
        [Writing a file called "fileName.txt" in saveLoc that
        contains initialString and then the timestamp.]
    '''

    from os.path import join

    outFile = join(saveLoc, fileName + '.txt')
    stringToWrite = initialString + now() + '\n'

    writeFile(outFile, stringToWrite, silent=True)
