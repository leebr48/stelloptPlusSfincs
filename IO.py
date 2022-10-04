# This file contains IO helper functions.

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
        [name1_AUX_S, name1_AUX_F, ...]. Note that the radial
        name (S) always comes before the value name (F).
    '''
  
    output_names = []
    for prefix in listOfPrefixes:
        s = prefix.lower().strip() + '_aux_s'
        f = prefix.lower() + '_aux_f'
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
        A list with two elements. The first is a list with all the unique 
        prefixes in nameList. The second is a list with interpolation 
        objects corresponding to the prefixes in the first list. The 
        independent (radial) variable for these interpolations is the 
        normalized toriodal flux. 
    '''
    
    import warnings
    
    strippedNames = []
    matched = []
    for namePair in nameList:
        strippedName = namePair[0].split('_')[0]
        strippedNames.append(strippedName)

        matchedPair = []
        for name in namePair:
            foundMatch = False
            for dataVec in dataList:
                if dataVec[0] == name:
                    matchedPair.append([float(i) for i in dataVec[1:]])
                    foundMatch = True
            if not foundMatch:
                warnings.warn('No match could be found for the variable "{}" in the given dataList!'.format(name))
        matched.append(matchedPair)

    if not any(matched):
        raise IOError('No searched variables were found.')
    
    return [strippedNames, matched]

def generatePreamble(radial_coordinate_ID):

    '''
    Inputs:
        The (integer) ID for the radial coordinate used in profiles.xxx. 0 = psiHat, 1 = psiN, 2 = rHat, 3 = rN.
    Outputs:
        String containing the preamble to the profiles.xxx file.
    '''

    stringToWrite = '# This is an integer specifying the radial coordinate used in this file, which can be different from the one specified by inputRadialCoordinate in input.namelist.\n'
    stringToWrite += '{}\n'.format(str(radial_coordinate_ID))
    stringToWrite += '# The following lines contain profile information in this format:\n'
    stringToWrite += '# radius\tNErs\tgeneralEr_min\tgeneralEr_max\tnHat(species 1)\tTHat(species 1)\tnHat(species 2)\tTHat(species 2)\t...\n'
    stringToWrite += '# The format of the generalEr_* profiles is set by inputRadialCoordinateForGradients in input.namelist. The default is Er.\n'
    stringToWrite += '\n'

    return stringToWrite
