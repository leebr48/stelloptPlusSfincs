# This file contains IO helper functions.

def listifyBEAMS3DFile(inputFile):
    
    '''
    Input:  STELLOPT input.namelist file with a BEAMS3D section.
    Output: The variable assignment data in the BEAMS3D section
            as a list of lists. Each sublist contains the variable
            name as its first element. Each subsequent element is 
            a value assigned to the variable. This means that 
            singly-valued variables will have lists with length 2
            and vector-valued variables will have lists with length 
            len(vector)+1.
    '''

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

    return listifiedData
