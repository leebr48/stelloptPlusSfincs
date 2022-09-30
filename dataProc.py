# This file contains functions for data processing.

def linearInterp(dataOfInterest):

    '''
    Inputs:
        dataOfInterest list, as from the extractDataList function.
    Outputs:
        The same list, but with SciPy interplation objects in place
        of the sublists in output[1].
    '''

    from scipy.interpolate import interp1d
    
    interps = []
    for dataPair in dataOfInterest[1]:
        interpObj = interp1d(dataPair[0],dataPair[1])
        interps.append(interpObj)

    return [dataOfInterest[0], interps]
