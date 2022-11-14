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
    parser.add_argument('--profilesIn', type=str, nargs='*', required=True, help='File(s) with relevant profiles written as in STELLOPT, with path(s) if necessary. This script currently reads the BEAMS3D section of STELLOPT namelist files. If you input multiple files, order matters!')
    parser.add_argument('--eqIn', type=str, nargs='*', help='File(s) from which to load the magnetic equilibrium(ia). Can be VMEC wout file(s) in netCDF or ASCII format, or IPP .bc file(s). If you input multiple files, order matters!')
    parser.add_argument("--vars", type=str, nargs='*', required=False, default=['NE', 'TI', 'NE', 'TE'], help='''Prefixes of variables to be read, normalized, and written. You should enter each prefix in quotes and put spaces between prefixes. The prefix names are not case sensitive. The density and temperature prefixes should come in the format <'N1' 'T1' 'N2' 'T2' ...> where '1' and '2' often indicate species identifiers (such as 'I' or 'E'). Note that you can write duplicate data by repeating entries. For instance, inputting <'NE' 'TI' 'NE' 'TE'> enforces NI=NE. The order in which the species prefixes are specified should match the species order in input.namelist. If you have potential data to input to calculate the radial electric field, 'POT' can be added anywhere in the list. The potential should give -Er when differentiated with respect to the STELLOPT coordinate S, which is psiN in SFINCS.''')
    parser.add_argument('--minBmn', type=float, nargs=1, required=False, default=[0.0], help='Only Fourier modes of at least this size will be loaded from the <eqIn> file(s).')
    parser.add_argument('--Nyquist', type=int, nargs=1, required=False, default=[2], help='Include the larger poloidal and toroidal mode numbers in the xm_nyq and xn_nyq arrays, where available, if this parameter is set to 2. Exclude these mode numbers if this parameter is set to 1.')
    parser.add_argument('--numInterpSurf', type=int, nargs=1, required=False, default=[1000], help='Number of radial surfaces on which to calculate and write interpolated profile data. This number should be quite large.')
    parser.add_argument('--radialVar', type=int, nargs=1, required=False, default=[3], help='ID of the radial coordinate used in the input.namelist file to specify which surfaces should be scanned over. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, and 3 = rN (which is the STELLOPT rho)')
    parser.add_argument('--radialGradientVar', type=int, nargs=1, required=False, default=[4], help='ID of the radial coordinate used to take derivatives. Relevant for the generalEr_* parameters in the profiles file and specifying the density and temperature derivatives on a single flux suface. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, 3 = rN (which is the STELLOPT rho), and 4 = rHat (like option 2, except that Er is used in place of dPhiHatdrHat). The default is recommended.')
    parser.add_argument('--numCalcSurf', type=int, nargs=1, required=False, default=[16], help='Number of radial surfaces on which to perform full SFINCS calculations.')
    parser.add_argument('--minRad', type=float, nargs=1, required=False, default=[0.15], help='Lower bound for the radial scan. If <resScan> is used, the flux surface specified by this parameter will be used for the convergence scan. Note that VMEC has resolution issues near the magnetic axis and SFINCS often converges much slower there due to the relatively low collisionality, so setting <minRad> to be very small may cause problems.')
    parser.add_argument('--maxRad', type=float, nargs=1, required=False, default=[0.95], help='Upper bound for the radial scan.')
    parser.add_argument('--Zs', type=float, nargs='*', required=False, default=[1, -1], help='Charge of each species in units of the proton charge. The species ordering must match that in the <vars> option.')
    parser.add_argument('--mHats', type=float, nargs='*', required=False, default=[1, 0.000545509], help='Mass of each species in units of the proton mass.')
    parser.add_argument('--seedEr', type=float, nargs=1, required=False, default=[0], help="Input an initial guess for the radial electric field in units of <radialGradientVar>. The default value is typically fine. This will be overwritten if you trigger an electric field scan with <numManErScan>.")
    parser.add_argument('--numManErScan', type=int, nargs=1, required=False, default=[0], help='Number of manual radial electric field scans to perform (using scanType=5). This parameter generates equidistant radial electric field seed values between <minEr> and <maxEr> for the root-finding algorithm in SFINCS. This parameter will be overwritten if <resScan> is activated.')
    parser.add_argument('--minEr', type=float, nargs=1, required=False, default=[-5], help='Minimum seed value of the radial electric field in units of <radialGradientVar>. This parameter is also used to derive the minimum and maximum Er available to ambipolarSolve. Note that you may need to change this to get good results. It is suggested that you seed ambipolarSolve with values near Er=0 initially, otherwise the solver could fail or converge to a very large (and erroneous) value of Er. Keep in mind that the value of the radial current at the "smallest electric field available to ambipolarSolve" (mentioned before) must have the opposite sign of the radial current at the "largest electric field available to ambipolarSolve" (mentioned in the <maxEr> help).')
    parser.add_argument('--maxEr', type=float, nargs=1, required=False, default=[5], help='Maximum seed value of the radial electric field in units of <radialGradientVar>. This parameter is also used to derive the minimum and maximum Er available to ambipolarSolve. Note that you may need to change this to get good results. It is suggested that you seed ambipolarSolve with values near Er=0 initially, otherwise the solver could fail or converge to a very large (and erroneous) value of Er. Keep in mind that the value of the radial current at the "largest electric field available to ambipolarSolve" (mentioned before) must have the opposite sign of the radial current at the "smallest electric field available to ambipolarSolve" (mentioned in the <minEr> help).')
    parser.add_argument('--resScan', action='store_true', default=False, help='Triggers a SFINCS resolution scan run.')
    parser.add_argument('--defaultDens', type=float, nargs='*', required=False, default=[1, 1], help='If <resScan> is used, this sets the density of each species in units of 1e20 m^-3. Note that you must specify a density for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultTemps', type=float, nargs='*', required=False, default=[1, 1], help='If <resScan> is used, this sets the temperature of each species in keV. Note that you must specify a temperature for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultDensDer', type=float, nargs='*', required=False, default=[-0.5e0, -0.5e0], help='If <resScan> is used, this sets the derivative of the density of each species (in units of 1e20 m^-3) with respect to psiN (which is the STELLOPT "S"). Note that you must specify a value for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultTempsDer', type=float, nargs='*', required=False, default=[-2e0, -2e0], help='If <resScan> is used, this sets the derivative of the temperature of each species (in keV) with respect to psiN (which is the STELLOPT "S"). Note that you must specify a value for EACH species. The exact values are probably not important.')
    parser.add_argument('--Nzeta', type=int, nargs=1, required=False, default=[55], help='Number of toroidal grid points per period. This should be an odd number.')
    parser.add_argument('--NzetaScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nzeta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Ntheta', type=int, nargs=1, required=False, default=[25], help='Number of poloidal grid points. This should be an odd number.')
    parser.add_argument('--NthetaScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Ntheta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nxi', type=int, nargs=1, required=False, default=[90], help='Number of Legendre polynomials used to represent the pitch-angle dependence of the distribution function.')
    parser.add_argument('--NxiScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nxi that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nx', type=int, nargs=1, required=False, default=[19], help='Number of grid points in energy used to represent the distribution function.')
    parser.add_argument('--NxScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nx that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--NL', type=int, nargs=1, required=False, default=[5], help='Number of Legendre polynomials used to represent the Rosenbluth potentials. Increasing this hardly changes the results, so it can almost certainly be left alone.')
    parser.add_argument('--NLScan', type=float, nargs=2, required=False, default=[0.5, 1.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of NL that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--solverTol', type=float, nargs=1, required=False, default=[1e-6], help='Tolerance used to define convergence of the iterative (Krylov) solver.')
    parser.add_argument('--solverTolScan', type=float, nargs=2, required=False, default=[0.1, 10.0], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of solverTolerance that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--saveLoc', type=str, nargs='*', required=False, default=[None], help='Location(s) in which to save written files - this will act as the main directory(ies) for a set of SFINCS runs. Defaults to either <profilesIn> or <eqIn> location(s) -- whichever input has more locations will be chosen as the default. If len(<profilesIn>) == len(<eqIn>), defaults to <profilesIn> location(s). If you input multiple files, order matters! Note that if you specify 1 <saveLoc> and multiple <profilesIn> or <eqIn>, the code will attempt to save all the generated files in the same directory. Due to the current (strict) naming conventions of SFINCS, this is probably not useful because the last-written files will overwrite their predecessors, but the feature is included for completeness.')
    parser.add_argument('--nNodes', type=int, nargs=1, required=False, default=[None], help='Total number of nodes to use for each SFINCS run. You must specify at least one of <nNodes> and <nTasks>.')
    parser.add_argument('--nTasksPerNode', type=int, nargs=1, required=False, default=[None], help='Number of MPI tasks to use on each node for each SFINCS run. This parameter should only be used if <nNodes> is specified and should not be used with <nTasks>.')
    parser.add_argument('--nTasks', type=int, nargs=1, required=False, default=[None], help='Total number of MPI tasks to use for each SFINCS run. You must specify at least one of <nNodes> and <nTasks>.')
    parser.add_argument('--mem', type=int, nargs=1, required=False, default=[None], help='Total amount of memory (MB) allocated for each SFINCS run.')
    parser.add_argument('--time', type=str, nargs=1, required=False, default=['00-18:00:00'], help='Wall clock time limit for the batch runs. Format is DD-HH:MM:SS. Note that SFINCS typically has the most trouble converging near the magnetic axis (due to the lower collisionality there cause by peaked temperature profiles), so you may need to increase <time> for runs near the axis.')
    parser.add_argument('--noProfiles', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write a profiles file.')
    parser.add_argument('--noNamelist', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write an input.namelist file.')
    parser.add_argument('--noBatch', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write a job.sfincsScan file.')
    parser.add_argument('--noRun', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not run sfincsScan.')
    parser.add_argument('--notifs', type=str, nargs=1, required=False, default=['bad'], help='Dictate which Slurm notification emails you would like to receive. By default, you will only receive emails when something bad happens to your job (such as a failure). You may also specify "all" or "none", which have the (intuitive) meanings indicated in the Slurm documentation. Note that the environment variable SFINCS_BATCH_EMAIL must be set for <notifs> to work correctly.')
    parser.add_argument('--noConfirm', action='store_true', default=False, help='Instruct sfincsScan to create folders and jobs without asking for confirmation first.')
    args = parser.parse_args()

    if args.minEr >= args.maxEr:
        raise IOError('<minEr> must be less than <maxEr>.')

    if args.Nyquist[0] not in [1,2]:
        raise IOError('An invalid <Nyquist> choice was specified. Valid inputs are the integers 1 and 2.')

    if args.radialVar[0] not in [0,1,2,3]:
        raise IOError('An invalid <radialVar> choice was specified. Valid inputs are the integers 0, 1, 2, and 3.')
    
    if args.radialGradientVar[0] not in [0,1,2,3,4]:
        raise IOError('An invalid <radialGradientVar> choice was specified. Valid inputs are the integers 0, 1, 2, 3, and 4.')

    if len(args.Zs) != len(args.vars)/2 and len(args.Zs) != (len(args.vars)-1)/2:
        raise IOError('The <Zs> input length is inconsistent with the <vars> input length.')
    
    if len(args.mHats) != len(args.vars)/2 and len(args.mHats) != (len(args.vars)-1)/2:
        raise IOError('The <mHats> input length is inconsistent with the <vars> input length.')
    
    if len(args.defaultDens) != len(args.vars)/2 and len(args.defaultDens) != (len(args.vars)-1)/2:
        raise IOError('The <defaultDens> input length is inconsistent with the <vars> input length.')
    
    if len(args.defaultTemps) != len(args.vars)/2 and len(args.defaultTemps) != (len(args.vars)-1)/2:
        raise IOError('The <defaultTemps> input length is inconsistent with the <vars> input length.')

    if len(args.defaultDensDer) != len(args.vars)/2 and len(args.defaultDensDer) != (len(args.vars)-1)/2:
        raise IOError('The <defaultDensDer> input length is inconsistent with the <vars> input length.')
    
    if len(args.defaultTempsDer) != len(args.vars)/2 and len(args.defaultTempsDer) != (len(args.vars)-1)/2:
        raise IOError('The <defaultTempsDer> input length is inconsistent with the <vars> input length.')

    if args.Nzeta[0]%2 == 0:
        raise IOError('<Nzeta> should be odd.')
    
    if args.Ntheta[0]%2 == 0:
        raise IOError('<Ntheta> should be odd.')
    
    lens = [len(args.profilesIn), len(args.eqIn), len(args.saveLoc)]
    maxLen = max(lens)
    for length in lens:
        if length != 1 and length != maxLen:
            raise IOError('Regarding <profilesIn>, <eqIn>, and <saveLoc>: any of these three inputs with length greater than 1 must have the same length as the other inputs with length greater than 1.')

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
    parser.add_argument('--sfincsDir', type=str, nargs='*', required=True, help='Top directory(ies) for SFINCS run(s), with path(s) if necessary. Such directories typically contain subdirectories which either contain SFINCS output files (*.h5) or more subdirectories for the electric field scan. In the latter case, those subsubdirectories contain SFINCS output files. If you input multiple directories, order matters!')
    parser.add_argument('--radialVar', type=int, nargs=1, required=False, default=[3], help='ID of the radial coordinate used in the input.namelist file to specify which surfaces should be scanned over. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, and 3 = rN (which is the STELLOPT rho)')
    parser.add_argument('--radialVarBounds', type=float, nargs=2, required=False, default=[-1, -1], help='Two floats, which are (in order) the minimum and maximum values of <radialVar> that will be plotted. If one of the inputs is negative, it will be ignored (so that the min or max is not limited).')
    parser.add_argument('--saveLoc', type=str, nargs='*', required=False, default=[None], help='Location(s) in which to save plots, plot data, and informational *.txt files. Defaults to <sfincsDir>/processed/. If you input multiple directories, order matters!')
    parser.add_argument('--checkConv', action='store_true', default=False, help='Instead of plotting anything, just check if the SFINCS runs in the <sfincsDir> location(s) converged. If they all did, you will receive no output.')
    args = parser.parse_args()

    if not all([isdir(item) for item in args.sfincsDir]):
        raise IOError('The inputs given in <sfincsDir> must be directories, not files.')

    if args.radialVar[0] not in [0,1,2,3]:
        raise IOError('An invalid <radialVar> choice was specified. Valid inputs are the integers 0, 1, 2, and 3.')
    
    lens = [len(args.sfincsDir), len(args.saveLoc)]
    maxLen = max(lens)
    for length in lens:
        if length != 1 and length != maxLen:
            raise IOError('If both <sfincsDir> and <saveLoc> have length greater than 1, they must be the same length.')
    
    return args

def getFileInfo(inFile, saveLoc, outFileName):

    '''
    Inputs:
        inFile: String with (relative or absolute) path
                to an input file.
        saveLoc: String with (relative or absolute) path 
                 where other files (such as the outputs
                 of other scripts) should be saved. Defaults
                 to the location of inFile, but can be
                 specified with <saveLoc> command line option.
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
    
    import itertools
        
    with open(inputFile,'r') as f:
        beams3dSectionStartFlag = False
        beams3dSectionEndFlag = False
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
            raise IOError('The script thinks that there were two equals signs in a variable assignment line of the input file. Something is wrong.')
        
        if (' ' or '\t') in precleaned[1]:
            cleaned = [precleaned[0]] + precleaned[1].split()
        else:
            cleaned = precleaned

        listifiedData.append(cleaned)

    listifiedData.sort()

    redundanciesRemoved = list(listifiedData for listifiedData,_ in itertools.groupby(listifiedData))

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
        A list of BEAMS3D variables using listOfPrefixes, such as
        [name1_AUX_S, name1_AUX_F, ...].
    '''
  
    output_names = []
    for prefix in listOfPrefixes:
        s = prefix + '_aux_s'
        f = prefix + '_aux_f'
        appendor = [s, f]
        output_names.append(appendor)

    return output_names

def extractDataList(dataList, nameList):

    '''
    Inputs:  
        dataList: A list of lists, as from the listifyBEAMS3DFile function,
                   with a string as the first element and length >= 2.
        nameList: A list of lists, as from the makeProfileNames function. 
                  Each sublist contains a pair of strings to search for in 
                  dataList.
    Outputs:
        A dictionary. Each key is a unique prefix from nameList. Each value
        contains a dictionary with two entries. The 'iv' entry is the radial
        coordinate (normalized toroidal flux) list for the given variable. 
        The 'dv' entry gives the list of values corresponding to those
        radial coordinates.
    '''

    import warnings
    
    matched = []
    dataDict = {}
    for namePair in nameList:
        
        strippedName = namePair[0].split('_')[0]

        matchedPair = {}
        for name in namePair:
            foundMatch = False

            for dataVec in dataList:
                
                if dataVec[0] == name:
                    floats = [float(i) for i in dataVec[1:]]
                    
                    if name[-1] == 's':
                        matchedPair['iv'] = floats
                    elif name[-1] == 'f':
                        matchedPair['dv'] = floats
                    else:
                        raise IOError('The read variable suffix is not "S" or "F". Something is wrong.')
                    
                    foundMatch = True

            if not foundMatch:
                warnings.warn('No match could be found for the variable "{}" in the given dataList!'.format(name))
        
        matched.append(matchedPair)
        dataDict[strippedName] = matchedPair

    if not any(matched):
        raise IOError('No searched variables were found.')
    
    return dataDict

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

def findFiles(name, path):

    '''
    Inputs:
        name: string with the name of a file to search for.
              Note that the name must be exact (there is
              no name matching).
        path: string with path of directory to search 
              (recursively) for name.
    Outputs:
        List with absolute paths to files called name
        within path.
    '''
    
    from os import walk
    from os.path import join

    result = []
    for root, dirs, files in walk(path):
        if name in files:
            result.append(join(root, name))
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
                outListDict
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
        
        extensiveParticleFluxUnits = r' $\mathrm{\left(\frac{1}{s}\right)}$'
        extensiveHeatFluxUnits = r' $\mathrm{\left(\frac{J}{s}\right)}$'
        extensiveMomentumFluxUnits = r' $\mathrm{\left(\frac{kg T m}{s^{-2}}\right)}$'
        extensiveRadialCurrentUnits = r' $\mathrm{\left(A\right)}$'
        
        if inString == 'Er':
            return r'Radial electric field $\mathrm{\left(\frac{V}{m}\right)}$'
        
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
            return r'Neoclassical heat flux' + extensiveHeatFluxUnits
        
        elif inString == 'extensiveClassicalHeatFlux':
            return r'Classical heat flux' + extensiveHeatFluxUnits
        
        elif inString == 'extensiveTotalHeatFlux':
            return r'Total heat flux' + extensiveHeatFluxUnits
        
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

        if len(parts) == 3 and parts[1] != 'vm': # With 'vm', the full distribution function is taken into account rather than only the leading-order contribution
            raise IOError('Formatting has not yet been specified for the variable {}.'.format(inString))
        
        label = parts[0]
        radVar = prettyRadialVar(parts[-1], innerOnly=True)
        
        # Write some label strings that will be used below
        directionStatement = r' in $\nabla {}$ direction '.format(radVar) # Note that this includes spaces on either side for convenience
        particleFluxUnits = r'$\mathrm{\left(\frac{1}{m^{3} s}\right)}$'
        heatFluxUnits = r'$\mathrm{\left(\frac{J}{m^{3} s}\right)}$'
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
            return r'Neoclassical heat flux' + directionStatement + heatFluxUnits
 
        elif label == 'classicalHeatFlux':
            return r'Classical heat flux' + directionStatement + heatFluxUnits
        
        elif label == 'classicalHeatFluxNoPhi1':
            return r'Classical heat flux (neglecting $\Phi_{1}$)' + directionStatement + heatFluxUnits

        elif label == 'totalHeatFlux':
            return r'Total heat flux' + directionStatement + heatFluxUnits

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
