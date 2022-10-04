# This file contains functions for data processing.

def scaleData(dataOfInterest, phiBar, nBar, TBar):

    '''
    Inputs:
        dataOfInterest: list, as from extractDataList function.
        phiBar: Reference value of the electrostatic potential in units of kV.
        nBar: Reference value of density in units of m^(-3).
        TBar: Reference value of temperature in units of keV.
    Outputs:
        dataOfInterest, but with the appropriate values scaled for SFINCS.
    '''
    
    import numpy as np

    scaledVals = []
    for i, name in enumerate(dataOfInterest[0]):
        
        radialData = dataOfInterest[1][i][0]
        variableData = dataOfInterest[1][i][1]
        
        if name[0] == 'p': # Potential
            multip = 1/1000/phiBar #BEAMS3D uses V instead of kV
        elif name[0] == 'n': # Density
            multip = 1/nBar
        elif name[0] == 't': # Temperature
            multip = 1/1000/TBar # BEAMS3D uses eV instead of keV
        else:
            raise IOError('I am not sure how to scale at least one of the input data arrays.')
            
        scaled = list(np.array(variableData)*multip)
        repaired = [radialData, scaled]
        scaledVals.append(repaired)

    return [dataOfInterest[0], scaledVals]

def linearInterp(inputData):

    '''
    Inputs:
        inputData list, as from the scaleData function.
    Outputs:
        inputData, but with SciPy interpolation objects in place
        of the sublists in output[1].
    '''

    from scipy.interpolate import interp1d
    
    interps = []
    for dataPair in inputData[1]:
        interpObj = interp1d(dataPair[0], dataPair[1])
        interps.append(interpObj)

    return [inputData[0], interps]

def nonlinearInterp(inputData, s=0, k=3, der=0):

    '''
    Inputs:
        inputData: list, as from the scaleData function.
        s: smoothing parameter. For details, see the docs on
           scipy.interpolate.splrep. s=0 corresponds to no
           smoothing, meaning that the interpolation function
           will go through every input data point.
        k: degree of the spline
        der: degree of derivative to be taken
    Outputs:
        inputData, but with SciPy interpolation objects in place
        of the sublists in output[1].
    '''

    from scipy.interpolate import splrep, splev

    interpObjects = []
    for dataPair in inputData[1]:
        tck = splrep(dataPair[0], dataPair[1], s=s, k=k)
        interpObj = lambda x, tck=tck, der=der: splev(x, tck, der, ext=2) # Will return error if extrapolation is requested
        interpObjects.append(interpObj)

    return [inputData[0], interpObjects]
