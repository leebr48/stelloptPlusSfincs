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

def findNumCalcs(baseVal, minMaxList, powersMode=False):

    '''
    Inputs:
        baseVal: base value of a calculation (such as a parameter
                 whose convergence is to be tested using sfincsScan).
        minMaxList: list containing the minimum and maximum values
                    by which baseVal is to be multiplied, such as in
                    sfincsScan.
        powersMode: set the output number of scans based on a logarithmic
                    scale rather than a linear scale. This is most 
                    appropriate for scanning over solver tolerances,
                    where the scan range will be very large.
    Outputs:
        A dictionary containing the minimum value in minMaxList (key
        'min'), the maximum value in minMaxList (key 'max'), and the
        (integer) maximum number of calculations that should be
        performed for the given variable during a sfincsScan operation
        (key 'num'). Note that this function currently gives a much 
        higher number of evaluations than will actually be used by 
        sfincsScan since many of the geometry parameters can only be
        odd integers.
    '''
    
    import numpy as np

    if len(minMaxList) != 2:
        raise IOError('An input list meant to contain minimum and maximum multipliers for a parameter to be scanned does not have two elements.')

    minMult = min(minMaxList)
    maxMult = max(minMaxList)

    if minMult == 0 and maxMult == 0:
        out = 0
    else:
        if powersMode:
            minimum = baseVal * minMult
            maximum = baseVal * maxMult
            out = int(np.ceil(np.log10(maximum / minimum)) + 1)
        else:
            out = int(np.ceil((maxMult - minMult) * baseVal))

    return {'min':minMult, 'max':maxMult, 'num':out}

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

    outputData = {}
    for key, data in inputData.items():
        tck = splrep(data['iv'], data['dv'], k=k, s=s)
        interpObj = lambda x, tck=tck, der=ders[key]: splev(x, tck, der=der, ext=2) # Will raise an error if extrapolation is requested
        outputData[key] = interpObj

    return outputData

def fixOutputUnits(inVar, inFloat, mBar=1.672621911e-27, BBar=1, RBar=1, nBar=1e20, TBar=1.60217733e-16, phiBar=1000):

    '''
    Inputs:
        inVar: name (string) of a SFINCS variable from its output (*.h5)
               file. Note that only a few variables are currently supported.
        inFloat: numerical value of inVar.
        mBar: reference mass. Default is the proton mass in kilograms.
              The default value is highly recommended.
        BBar: reference magnetic field in tesla.
              The default value is highly recommended.
        RBar: reference length in meters.
              The default value is highly recommended.
        nBar: reference density in meters^-3.
              The default value is highly recommended.
        TBar: reference temperature in joules. Default is 1 keV.
              The default value is highly recommended.
        phiBar: reference electric potential in volts. Default
                is 1 kV. The default value is highly recommended.
    Outputs:
        The value of inVar in SI units. 
    '''
    
    import numpy as np

    # Other important quantities
    e = 1.602176634e-19 # C (proton charge)
    vBar = np.sqrt(2 * TBar / mBar) # meters/second by default

    # Identify inVar and perform some administrative checks
    if '_' not in inVar: # These are not radial fluxes
        shouldHaveUnits = inVar
    else: # These are radial fluxes
        shouldHaveUnits = inVar.split('_')[0]

    # Convert to SI units 
    if shouldHaveUnits == 'Er':
        return phiBar / RBar * inFloat # V/m
    
    elif shouldHaveUnits == 'FSABjHat':
        return e * nBar * vBar * inFloat # A*m^-2, note that we do not multiply by BBar
    
    elif shouldHaveUnits == 'FSABFlow':
        return nBar * vBar * inFloat # m^-2*s^-1, note that we do not multiply by BBar
    
    elif shouldHaveUnits == 'particleFlux': # Neoclassical
        return nBar * vBar * inFloat # m^-2*s^-1, note that we do not divide by RBar
    
    elif shouldHaveUnits == 'extensiveParticleFlux': # Neoclassical
        return nBar * vBar * RBar**2 * inFloat # s^-1
    
    elif shouldHaveUnits in ['classicalParticleFlux', 'classicalParticleFluxNoPhi1', 'totalParticleFlux']:
        return nBar * vBar * inFloat # m^-2*s^-1, note that we do not multiply by Z*e or divide by RBar

    elif shouldHaveUnits == 'heatFlux': # Neoclassical
        return mBar * nBar * vBar**3 * inFloat # J*m^-2*s^-1 = kg*s^-3, note that we do not divide by RBar

    elif shouldHaveUnits == 'extensiveHeatFlux': # Neoclassical
        return mBar * nBar * vBar**3 * RBar**2 * inFloat # J*s^-1
    
    elif shouldHaveUnits in ['classicalHeatFlux', 'classicalHeatFluxNoPhi1', 'totalHeatFlux']:
        return mBar * nBar * vBar**3 * inFloat # J*m^-2*s^-1 = kg*s^-3, note that we do not multiply by Z*e or divide by RBar
    
    elif shouldHaveUnits == 'momentumFlux': # Neoclassical 
        return mBar * nBar * vBar**2 * inFloat # kg*m^-1*s^-2, note that we do not divide by RBar or multiply by BBar

    elif shouldHaveUnits == 'radialCurrent':
        return e * nBar * vBar * inFloat # A*m^-2

    else:
        raise IOError('Conversion factor has not yet been specified for the variable {}.'.format(inVar))
