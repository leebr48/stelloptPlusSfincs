# This file contains IO helper functions.

def getArgs():

    '''
    Inputs:
        [No direct inputs. See below for command line inputs.]
    Outputs:
        Arguments that can be passed to other scripts.
    '''

    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('profilesIn', type=str, nargs=1, help='File with relevant profiles written as in STELLOPT, with path if necessary. This script currently reads the BEAMS3D section of the STELLOPT namelist file.')
    parser.add_argument('eqIn', type=str, nargs=1, help='File from which to load the magnetic equilibrium. Can be a VMEC wout file in netCDF or ASCII format, or an IPP .bc file.')
    parser.add_argument("--vars", type=str, nargs='*', required=False, default=['NE', 'TI', 'NE', 'TE'], help='''Prefixes of variables to be read, normalized, and written. You should enter each prefix in quotes and put spaces between prefixes. The prefix names are not case sensitive. The density and temperature prefixes should come in the format <'N1' 'T1' 'N2' 'T2' ...> where '1' and '2' often indicate species identifiers (such as 'I' or 'E'). Note that you can write duplicate data by repeating entries. For instance, inputting <'NE' 'TI' 'NE' 'TE'> enforces NI=NE. The order in which the species prefixes are specified should match the species order in input.namelist. If you have potential data to input to calculate the radial electric field, 'POT' can be added anywhere in the list. The potential should give -Er when differentiated with respect to the STELLOPT coordinate S, which is psiN in SFINCS.''')
    parser.add_argument('--minBmn', type=float, nargs=1, required=False, default=[0.0], help='Only Fourier modes of at least this size will be loaded from the <eqIn> file.')
    parser.add_argument('--Nyquist', type=int, nargs=1, required=False, default=[2], help='Include the larger poloidal and toroidal mode numbers in the xm_nyq and xn_nyq arrays, where available, if this parameter is set to 2. Exclude these mode numbers if this parameter is set to 1.')
    parser.add_argument('--numInterpSurf', type=int, nargs=1, required=False, default=[1000], help='Number of radial surfaces on which to calculate and write interpolated profile data. This number should be quite large.')
    parser.add_argument('--radialVar', type=int, nargs=1, required=False, default=[3], help='ID of the radial coordinate used in the input.namelist file to specify which surfaces should be scanned over. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, and 3 = rN (which is the STELLOPT rho)')
    parser.add_argument('--radialGradientVar', type=int, nargs=1, required=False, default=[4], help='ID of the radial coordinate used to take derivatives. Relevant for the generalEr_* parameters in the profiles file and specifying the density and temperature derivatives on a single flux suface. Valid entries are: 0 = psiHat, 1 = psiN (which is the STELLOPT "S"), 2 = rHat, 3 = rN (which is the STELLOPT rho), and 4 = rHat (like option 2, except that Er is used in place of dPhiHatdrHat). The default is recommended.')
    parser.add_argument('--numCalcSurf', type=int, nargs=1, required=False, default=[16], help='Number of radial surfaces on which to perform full SFINCS calculations.')
    parser.add_argument('--minRad', type=float, nargs=1, required=False, default=[0.05], help='Lower bound for the radial scan. If <resScan> is used, the flux surface specified by this parameter will be used for the convergence scan.')
    parser.add_argument('--maxRad', type=float, nargs=1, required=False, default=[0.95], help='Upper bound for the radial scan.')
    parser.add_argument('--Zs', type=float, nargs='*', required=False, default=[1, -1], help='Charge of each species in units of the proton charge. The species ordering must match that in the <vars> option.')
    parser.add_argument('--mHats', type=float, nargs='*', required=False, default=[1, 0.000545509], help='Mass of each species in units of the proton mass.')
    parser.add_argument('--seedEr', type=float, nargs=1, required=False, default=[0], help="Input an initial guess for the radial electric field in units of <radialGradientVar>. The default value is typically fine. This will be overwritten if you trigger an electric field scan with <numManErScan>.")
    parser.add_argument('--numManErScan', type=int, nargs=1, required=False, default=[0], help='Number of manual radial electric field scans to perform (using scanType=5). Note that the internal root-finding algorithms in SFINCS are quite reliable and will usually converge to the correct answer. If they fail, you can give them better initial guesses using this parameter, <minEr>, and <maxEr>. This parameter will be overwritten if <resScan> is activated.')
    parser.add_argument('--minEr', type=float, nargs=1, required=False, default=[-50], help='If a radial electric field scan should occur: minimum seed value of the generalized Er variable. Note that you may need to change this to get good results.')
    parser.add_argument('--maxEr', type=float, nargs=1, required=False, default=[50], help='If a radial electric field scan should occur: maximum seed value of the generalized Er variable. Note that you may need to change this to get good results.')
    parser.add_argument('--resScan', action='store_true', default=False, help='Triggers a SFINCS resolution scan run.')
    parser.add_argument('--defaultDens', type=float, nargs='*', required=False, default=[1, 1], help='If <resScan> is used, this sets the density of each species in units of 1e20 m^-3. Note that you must specify a density for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultTemps', type=float, nargs='*', required=False, default=[1, 1], help='If <resScan> is used, this sets the temperature of each species in keV. Note that you must specify a temperature for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultDensDer', type=float, nargs='*', required=False, default=[-0.5e0, -0.5e0], help='If <resScan> is used, this sets the derivative of the density of each species (in units of 1e20 m^-3) with respect to psiN (which is the STELLOPT "S"). Note that you must specify a value for EACH species. The exact values are probably not important.')
    parser.add_argument('--defaultTempsDer', type=float, nargs='*', required=False, default=[-2e0, -2e0], help='If <resScan> is used, this sets the derivative of the temperature of each species (in keV) with respect to psiN (which is the STELLOPT "S"). Note that you must specify a value for EACH species. The exact values are probably not important.')
    parser.add_argument('--Nzeta', type=int, nargs=1, required=False, default=[45], help='Number of toroidal grid points per period. This should be an odd number.')
    parser.add_argument('--NzetaScan', type=float, nargs=2, required=False, default=[0.5, 2.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nzeta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Ntheta', type=int, nargs=1, required=False, default=[45], help='Number of poloidal grid points. This should be an odd number.')
    parser.add_argument('--NthetaScan', type=float, nargs=2, required=False, default=[0.5, 2.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Ntheta that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nxi', type=int, nargs=1, required=False, default=[75], help='Number of Legendre polynomials used to represent the pitch-angle dependence of the distribution function.')
    parser.add_argument('--NxiScan', type=float, nargs=2, required=False, default=[0.5, 2.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nxi that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--Nx', type=int, nargs=1, required=False, default=[9], help='Number of grid points in energy used to represent the distribution function.')
    parser.add_argument('--NxScan', type=float, nargs=2, required=False, default=[0.5, 2.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of Nx that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--NL', type=int, nargs=1, required=False, default=[5], help='Number of Legendre polynomials used to represent the Rosenbluth potentials. Increasing this hardly changes the results, so it can almost certainly be left alone.')
    parser.add_argument('--NLScan', type=float, nargs=2, required=False, default=[0.5, 2.5], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of NL that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--solverTol', type=float, nargs=1, required=False, default=[1e-6], help='Tolerance used to define convergence of the iterative (Krylov) solver.')
    parser.add_argument('--solverTolScan', type=float, nargs=2, required=False, default=[0.1, 10.0], help='Two floats, which are (in order) the minimum and maximum multipliers on the value of solverTolerance that will be used if a resolution scan is run. Set both values to zero to not scan this parameter.')
    parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=[None], help='Location in which to save written files - this will act as the main directory for a set of SFINCS runs. Defaults to <profilesIn> location.')
    parser.add_argument('--nTasks', type=int, nargs=1, required=False, default=[8], help='Total number of MPI tasks to use for the SFINCS runs.')
    parser.add_argument('--mem', type=int, nargs=1, required=False, default=[75000], help='Total amount of memory (MB) allocated for the SFINCS runs.')
    parser.add_argument('--time', type=str, nargs=1, required=False, default=['00-00:10:00'], help='Wall clock time limit for the batch runs. Format is DD-HH:MM:SS.')
    parser.add_argument('--noProfiles', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write a profiles file.')
    parser.add_argument('--noNamelist', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write an input.namelist file.')
    parser.add_argument('--noBatch', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not write a job.sfincsScan file.')
    parser.add_argument('--noRun', action='store_true', default=False, help='Instruct higher-level wrapper scripts to not run sfincsScan.')
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

    import os

    inFile = os.path.abspath(inFile)
    inFilePath = os.path.dirname(inFile)
    inFileName = os.path.basename(inFile)

    if saveLoc == None:
        outFilePath = inFilePath
    else:
        outFilePath = os.path.abspath(saveLoc)

    outFile = os.path.join(outFilePath, outFileName)

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

def writeFile(outFile, stringToWrite):

    '''
    Inputs:
        outFile: absolute name (with path) of file
                 to be written.
        stringToWrite: string to write in outFile.
    Outputs:
        [A written file and a notification message.]
    '''
   
    _, outFileName, _, _, _ = getFileInfo(outFile, 'arbitrary/path/', 'arbitrary')

    with open(outFile, 'w') as f:
        f.write(stringToWrite)

    print('{} file written.'.format(outFileName))
