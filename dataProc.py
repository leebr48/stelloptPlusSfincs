# This file contains functions for data processing.

def findMinMax(dataOfInterest):

    '''
    Inputs:
        Dictionary, as from the extractDataList function.
    Outputs:
        Dictionary that indicates the largest minimum bound
        and smallest maximum bound from the independent
        variables in dataOfInterest. This information 
        would typically be used to determine the appropriate
        radial domain over which all the profile information
        could be accurately specified.
    '''
    
    out = {}
    mins = []
    maxes = []
    for key,data in dataOfInterest.items():
        mins.append(min(data['iv']))
        maxes.append(max(data['iv']))

    out['min'] = max(mins)
    out['max'] = min(maxes)

    return out

def scaleData(dataOfInterest, phiBar, nBar, TBar):

    '''
    Inputs:
        dataOfInterest: dictionary, as from the extractDataList function.
        phiBar: Reference value of the electrostatic potential in units of kV.
        nBar: Reference value of density in units of m^(-3).
        TBar: Reference value of temperature in units of keV.
    Outputs:
        dataOfInterest, but with the appropriate values scaled for SFINCS.
    '''
    
    import numpy as np

    for key, data in dataOfInterest.items():
        
        if key[0] == 'p': # Potential
            multip = 1/1000/phiBar #BEAMS3D uses V instead of kV
        elif key[0] == 'n': # Density
            multip = 1/nBar
        elif key[0] == 't': # Temperature
            multip = 1/1000/TBar # BEAMS3D uses eV instead of keV
        else:
            raise IOError('I am not sure how to scale at least one of the input data arrays.')
            
        scaled = list(np.array(data['dv'])*multip)

        dataOfInterest[key]['dv'] = scaled

    return dataOfInterest

def nonlinearInterp(inputData, s=0, k=3, der=0):

    '''
    Inputs:
        inputData: dictionary, as from the scaleData function.
        s: smoothing parameter. For details, see the docs on
           scipy.interpolate.splrep. s=0 corresponds to no
           smoothing, meaning that the interpolation function
           will go through every data point.
        k: degree of the spline
        der: order of derivative to be taken
    Outputs:
        inputData, but with SciPy interpolation objects in place
        of the sublists in output[1].
    '''

    from scipy.interpolate import splrep, splev

    for key, data in inputData.items():
        tck = splrep(data['iv'], data['dv'], s=s, k=k)
        interpObj = lambda x, tck=tck, der=der: splev(x, tck, der=der, ext=2) # Will return error if extrapolation is requested
        inputData[key] = interpObj

    return inputData
