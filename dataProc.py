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

def scaleData(dataOfInterest, phiBar=1, nBar=1e20, TBar=1):

    '''
    Inputs:
        dataOfInterest: Dictionary, as from the extractDataList function.
        phiBar: Reference value of the electrostatic potential in units of kV.
                Changing this will break other things!
        nBar: Reference value of density in units of m^(-3). Note that Python
              "e" notation is equivalent to Fortran "d" notation.
              Changing this will break other things!
        TBar: Reference value of temperature in units of keV.
              Changing this will break other things!
    Outputs:
        dataOfInterest, but with the appropriate values scaled for SFINCS.
    '''

    for key, data in dataOfInterest.items():
        
        if key[0] == 'p': # Potential
            multip = 1/1000/phiBar #BEAMS3D uses V instead of kV
        elif key[0] == 'n': # Density
            multip = 1/nBar
        elif key[0] == 't': # Temperature
            multip = 1/1000/TBar # BEAMS3D uses eV instead of keV
        else:
            raise IOError('I am not sure how to scale at least one of the input data arrays.')

        scaled = [item * multip for item in data['dv']]

        dataOfInterest[key]['dv'] = scaled

    return dataOfInterest

def nonlinearInterp(inputData, ders, k=3, s=0):

    '''
    Inputs:
        inputData: Dictionary, as from the scaleData function.
        ders: Dictionary with the same keys as inputData and
              values describing how many derivatives should
              be taken to determine the final interpolation
              object. The entries are typically zero, but may
              be 1 for calculating (say) derivatives of the 
              electric potential.
        k: Degree of the spline (prior to any derivatives 
           being taken).
        s: Smoothing parameter. For details, see the docs on
           scipy.interpolate.splrep. s=0 corresponds to no
           smoothing, meaning that the interpolation function
           will go through every data point.
    Outputs:
        inputData, but with SciPy interpolation objects in place
        of the data.
    '''

    from scipy.interpolate import splrep, splev

    for key, data in inputData.items():
        tck = splrep(data['iv'], data['dv'], k=k, s=s)
        interpObj = lambda x, tck=tck, der=ders[key]: splev(x, tck, der=der, ext=2) # Will return error if extrapolation is requested
        inputData[key] = interpObj

    return inputData
