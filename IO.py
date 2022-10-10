# This file contains IO helper functions.

def getArgs():

    '''
    Inputs:
        No direct inputs. See below for command line inputs.
    Outputs:
        Arguments that can be passed to other scripts.
    '''

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', type=str, nargs=1, required=True, help='File with relevant profiles written as in STELLOPT, with path if necessary. This script currently reads the BEAMS3D section of the STELLOPT namelist file.')
    parser.add_argument("--vars", type=str, nargs='*', required=True, help='''Prefixes of variables to be read, normalized, and written. You should enter each prefix in quotes and put spaces between prefixes. The prefix names are not case sensitive. The density and temperature prefixes should come in the format <'N1' 'T1' 'N2' 'T2' ...> where '1' and '2' often indicate species identifiers (such as 'I' or 'E'). Note that you can write duplicate data by repeating entries. For instance, inputting <'NE' 'TI' 'NE' 'TE'> enforces NI=NE. The order in which the species prefixes are specified should match the species order in input.namelist. If you have potential data to input to calculate the radial electric field, 'POT' can be added anywhere in the list. The potential should give -Er when differentiated with respect to the STELLOPT coordinate S, which is psiN in SFINCS.''')
    parser.add_argument('--saveLoc', type=str, nargs=1, required=False, default=None, help='Location in which to save profiles. Defaults to <inFile> location.')
    parser.add_argument('--numRad', type=int, nargs=1, required=False, default=[1000], help='Number of radial surfaces on which to calculate and write interpolated profile data. This number should be quite large. Default = 1000.')
    parser.add_argument('--numErScan', type=int, nargs=1, required=False, default=[5], help='If a radial electric field scan should occur: number of scans to perform. This parameter will be overwritten if Er data is provided. Default = 5.')
    parser.add_argument('--minEr', type=float, nargs=1, required=False, default=[-10], help='If a radial electric field scan should occur: minimum value of the generalized Er variable. Note that you may need to change this to get good results. This parameter will be overwritten if Er data is provided. Default = -10.')
    parser.add_argument('--maxEr', type=float, nargs=1, required=False, default=[10], help='If a radial electric field scan should occur: maximum value of the generalized Er variable. Note that you may need to change this to get good results. This parameter will be overwritten if Er data is provided. Default = 10.')
    parser.add_argument('--constEr', action='store_true', required=False, help='Assume the radial electric field is constant (as in scanType = 4). Er is set in input.namelist in this case.')
    parser.add_argument('--phiBar', type=float, nargs=1, required=False, default=[1], help='Reference electrostatic potential in units of kV. Default = 1.')
    parser.add_argument('--nBar', type=float, nargs=1, required=False, default=[1e20], help='Reference density in units of m^(-3). Note that Python "E" notation is equivalent to Fortran "D" notation. Default = 1e20.')
    parser.add_argument('--TBar', type=float, nargs=1, required=False, default=[1], help='Reference temperature in units of keV. Default = 1.')
    args = parser.parse_args()

    return args

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

    stringToWrite = '# This is an integer specifying the radial coordinate used in this file, which can be different from the one specified by inputRadialCoordinate in input.namelist.\n'
    stringToWrite += '{}\n'.format(str(radial_coordinate_ID))
    stringToWrite += '# The following lines contain profile information in this format:\n'
    stringToWrite += '# radius\tNErs\tgeneralEr_min\tgeneralEr_max\tnHat(species 1)\tTHat(species 1)\tnHat(species 2)\tTHat(species 2)\t...\n'
    stringToWrite += '# The meaning of the generalEr_* variables is set by inputRadialCoordinateForGradients in input.namelist. The default is Er_*.\n'

    return stringToWrite

def generateDataText(radii, *funcs):

    '''
    Inputs:
        radii: A list of radii at which to evaulate *funcs.
        *funcs: functions that can take one element of radii
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
