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
    Outputs:
        A list of BEAMS3D variables using listOfPrefixes, such as
        [name1_AUX_S, name1_AUX_F, ...]. Note that the radial
        name (S) always comes before the value name (F).
    '''
  
    import itertools

    redundanciesRemoved = list(listOfPrefixes for listOfPrefixes,_ in itertools.groupby(listOfPrefixes))

    output_names = []
    for prefix in redundanciesRemoved:
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
        A list with two elements. The first is nameList. The second is a 
        list of lists of lists. Within this second item, each sublist 
        contains a pair of lists that correspond to a pair of strings 
        passed to this function in nameList. This function preserves
        order - for instance, the name in output[0][0][0] corresponds
        to the data in output[1][0][0].
    '''
    
    import warnings
    
    matched = []
    for namePair in nameList:
        matchedPair = []
        for name in namePair:
            foundMatch = False
            for dataVec in dataList:
                if dataVec[0] == name:
                    matchedPair.append(dataVec[1:])
                    foundMatch = True
            if not foundMatch:
                warnings.warn('No match could be found for the variable "{}" in the given dataList!'.format(name))
        matched.append(matchedPair)

    if not any(matched):
        raise IOError('No searched variables were found.')
    
    return [nameList,matched]
