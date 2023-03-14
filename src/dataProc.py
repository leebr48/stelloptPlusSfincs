# This file contains functions for data processing.

def findMinMax(dataOfInterest):

    '''
    Inputs:
        Dictionary, as from the extractProfileData function.
    Outputs:
        Dictionary that indicates the largest minimum bound
        and smallest maximum bound from the independent
        variables in dataOfInterest. This information 
        would typically be used to determine the appropriate
        radial domain over which all the profile information
        could be accurately specified.
    '''

    import numpy as np
    
    out = {}
    mins = []
    maxes = []
    for key,data in dataOfInterest.items():
        mins.append(np.min(data['iv']))
        maxes.append(np.max(data['iv']))

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

def scaleInputData(dataOfInterest, profiles=True, phiBar=1, nBar=1e20, TBar=1, mBar=1.67262192369e-27, ZBar=1):

    '''
    Inputs:
        dataOfInterest: Dictionary, as from the extractProfileData or
                        extractScalarData functions.
        profiles: Boolean toggle. If True, dataOfInterest is assumed to contain
                  profile information (such as from extractProfileData). If
                  False, dataOfInterest is assumed to contain scalar information
                  (such as from extractScalarData).
        phiBar: Reference value of the electrostatic potential in units of kV.
                Changing this will break other things!
        nBar: Reference value of density in units of m^(-3). Note that Python
              "e" notation is equivalent to Fortran "d" notation.
              Changing this will break other things!
        TBar: Reference value of temperature in units of keV.
              Changing this will break other things!
        mBar: Reference value of mass (which is the proton mass) in kg.
              Changing this will break other things!
        ZBar: Reference value of charge in units of proton charge.
              Changing this will break other things!
    Outputs:
        dataOfInterest, but with the appropriate values scaled for SFINCS input.
    '''

    import numpy as np

    for key, data in dataOfInterest.items():
        
        if key[0] == 'p': # Potential
            multip = 1/1000/phiBar #BEAMS3D uses V instead of kV
        elif key[0] == 'n': # Density
            multip = 1/nBar
        elif key[0] == 't': # Temperature
            multip = 1/1000/TBar # BEAMS3D uses eV instead of keV
        elif key[0] == 'm': # Mass
            multip = 1/mBar # BEAMS3D uses kg instead of proton mass
        elif key[0] == 'z': # Charge
            multip = 1
        else:
            raise IOError('I am not sure how to scale at least one of the input data arrays.')
        
        if profiles:
            dataOfInterest[key]['dv'] = (multip * np.array(data['dv'])).tolist()
        else:
            dataOfInterest[key] = (multip * np.array(data)).tolist()

    return dataOfInterest

def constructBSpline(ivVec, dvVec, k=3, s=0):

    '''
    Inputs:
        ivVec: List of independent variable values.
        dvVec: List of dependent variable values.
        k: Degree of the spline.
        s: Smoothing parameter. For details, see the docs on
           scipy.interpolate.splrep. s=0 corresponds to no
           smoothing, meaning that the interpolation function
           will go through every data point.
    Outputs:
        A "tck" tuple describing the spline - see the docs
        mentioned above for details.
    '''
    
    from scipy.interpolate import splrep

    tck = splrep(ivVec, dvVec, k=k, s=s)
    
    return tck

def nonlinearInterp(inputData, ders, k=3, s=0, pchip=False):

    '''
    Inputs:
        inputData: Dictionary, as from the scaleInputData
                   function.
        ders: Dictionary with the same keys as inputData and
              values describing how many derivatives should
              be taken to determine the final interpolation
              object. The entries are typically zero, but may
              be 1 for calculating (say) derivatives of the 
              electric potential.
        k: Degree of the spline.
        s: Smoothing parameter. For details, see the docs on
           scipy.interpolate.splrep. s=0 corresponds to no
           smoothing, meaning that the interpolation function
           will go through every data point.
        pchip: If True, the SciPy PchipInterpolator will
               be used rather than a BSpline.
    Outputs:
        inputData, but with lists of SciPy interpolation
        objects in place of the data.
    '''

    from scipy.interpolate import splev
    from scipy.interpolate import PchipInterpolator

    outputData = {}
    for key, data in inputData.items():
        
        interpObjs = []
        for ivVec,dvVec in zip(data['iv'], data['dv']):
            if not pchip:
                tck = constructBSpline(ivVec, dvVec, k=k, s=s)
                interpObj = lambda x, tck=tck, der=ders[key]: splev(x, tck, der=der, ext=2) # Will raise an error if extrapolation is requested
            else:
                f = PchipInterpolator(ivVec, dvVec, extrapolate=False) # Will raise an error if extrapolation is requested
                interpObj = f.derivative(nu=ders[key])
            interpObjs.append(interpObj)
        
        outputData[key] = interpObjs

    return outputData

def fixOutputUnits(inVar, inFloat, mBar=1.67262192369e-27, BBar=1, RBar=1, nBar=1e20, TBar=1.602176634e-16, phiBar=1000):

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

    from scipy.constants import e # elementary charge in Coulombs

    # Important quantity
    vBar = thermalVelocity(TBar, mBar) # m/s by default

    # Identify inVar
    if '_' not in inVar:
        shouldHaveUnits = inVar
    else:
        shouldHaveUnits = inVar.split('_')[0]

    # Convert to SI units 
    # In the intensive radial fluxes, you will notice what appears to be an extra factor of m^-1.
    # This is not an error. It comes about because the heat and particle flux vectors carry m^-2
    # in their units, but they are always dotted with the gradient of a normalized radial
    # coordinate. The gradient itself carries a m^-1 unit, and because the radial coordinates are
    # normalized, they do not have any units. For presentations, it is suggested that you either
    # use extensive units for fluxes or multiply the intensive values by RBar to recover the
    # standard m^-2 units.
    
    if shouldHaveUnits == 'Er':
        return phiBar / RBar * inFloat # V/m
    
    elif shouldHaveUnits == 'FSABHat2':
        return BBar * inFloat # T^2

    elif shouldHaveUnits == 'FSABFlow':
        return nBar * vBar * BBar * inFloat # T*m^-2*s^-1
    
    elif shouldHaveUnits == 'FSABjHat':
        return e * nBar * vBar * BBar * inFloat # T*A*m^-2

    elif shouldHaveUnits in ['FSABjHatOverRootFSAB2', 'FSABjHatOverB0']:
        return e * nBar * vBar * inFloat # A*m^-2
     
    elif shouldHaveUnits in ['particleFlux', 'classicalParticleFlux', 'classicalParticleFluxNoPhi1', 'totalParticleFlux']:
        return nBar * vBar / RBar * inFloat # m^-3*s^-1
    
    elif shouldHaveUnits in ['extensiveParticleFlux', 'extensiveClassicalParticleFlux', 'extensiveTotalParticleFlux']:
        return nBar * vBar * RBar**2 * inFloat # s^-1
    
    elif shouldHaveUnits in ['heatFlux', 'classicalHeatFlux', 'classicalHeatFluxNoPhi1', 'totalHeatFlux']:
        return mBar * nBar * vBar**3 / RBar * inFloat # J*m^-3*s^-1

    elif shouldHaveUnits in ['extensiveHeatFlux', 'extensiveClassicalHeatFlux', 'extensiveTotalHeatFlux']:
        return mBar * nBar * vBar**3 * RBar**2 * inFloat # J*s^-1
    
    elif shouldHaveUnits == 'momentumFlux': # Neoclassical 
        return mBar * nBar * vBar**2 * BBar / RBar * inFloat # kg*T*m^-2*s^-2

    elif shouldHaveUnits == 'extensiveMomentumFlux': # Neoclassical
        return mBar * nBar * vBar**2 * BBar * RBar**2 * inFloat # kg*T*m*s^-2

    elif shouldHaveUnits == 'radialCurrent':
        return e * nBar * vBar / RBar * inFloat # A*m^-3

    elif shouldHaveUnits == 'extensiveRadialCurrent':
        return e * nBar * vBar * RBar**2 * inFloat # A

    else:
        raise IOError('Conversion factor has not yet been specified for the variable {}.'.format(inVar))

def convertRadDer(inputDerID, inputDerVal, outputDerID, aHat, psiAHat, psiN, XisPhi=False):

    '''
    Inputs:
        inputDerID: Integer specifying the radial variable with respect to which a derivative
                    is being taken. These values are specified in <radialGradientVar>.
        inputDerVal: Float specifying the value of the input derivative.
        outputDerID: Integer specifying the desired radial variable with respect to which
                     a derivative should be taken. These values are specified in
                     <radialGradientVar>.
        aHat: Float specifying the normalized effective minor radius at the last closed
              flux surface.
        psiAHat: Float specifying the normalized toroidal flux at the last closed flux
                 surface divided by 2*pi.
        psiN: Float specifying the toroidal flux normalized by its value at the last closed
              flux surface; equivalent to STELLOPT "s".
        XisPhi: Boolean specifying if the function that is being differentiated is the
                electric potential (True) or not (False); this is relevant when
                inputDerID or outputDerID are 4, because the radial electric field is used
                directly in this case rather than a derivative of the potential (which
                amounts to a change in the sign of the function output).
    Outputs:
        A float that is inputDerVal converted such that the derivative is taken
        with respect to outputDerID.
    '''

    from math import sqrt

    # These will save a bit of repetitive code below
    conv1 = 1 / psiAHat
    conv3 = conv1 / 2 / sqrt(psiN)
    conv2Or4 = conv3 * aHat

    # First convert input to dX/dpsiHat
    if inputDerID == 0: # dX/dpsiHat
        dXdpsiHat = inputDerVal
    elif inputDerID == 1: # dX/dpsiN
        dXdpsiHat = inputDerVal * conv1
    elif inputDerID == 2: # dX/drHat
        dXdpsiHat = inputDerVal * conv2Or4
    elif inputDerID == 3: # dX/drN
        dXpsiHat = inputDerVal * conv3
    elif inputDerID == 4 and not XisPhi: #dX/drHat
        dXdpsiHat = inputDerVal * conv2Or4
    elif inputDerID == 4 and XisPhi: #dX/drHat, except Er = -dPhiHat/drHat is used instead of dPhiHat/drHat
        dXdpsiHat = -1 * inputDerVal * conv2Or4
    else:
        raise IOError('An unknown inputDerID was passed to this function.')

    # Now convert dX/dpsiHat to the desired output
    if outputDerID == 0: # dX/dpsiHat
        outputDerVal = dXdpsiHat
    elif outputDerID == 1: # dX/dpsiN
        outputDerVal = dXdpsiHat / conv1
    elif outputDerID == 2: # dX/drHat
        outputDerVal = dXdpsiHat / conv2Or4
    elif outputDerID == 3: #dX/drN
        outputDerVal = dXdpsiHat / conv3
    elif outputDerID == 4 and not XisPhi: #dX/drHat
        outputDerVal = dXdpsiHat / conv2Or4
    elif outputDerID == 4 and XisPhi: #dX/drHat, except Er = -dPhiHat/drHat is used instead of dPhiHat/drHat
        outputDerVal = -1 * dXdpsiHat / conv2Or4
    else:
        raise IOError('An unknown outputDerID was passed to this function.')

    return outputDerVal

def checkConvergence(file):

    '''
    Inputs:
       Absolute path to a SFINCS output (*.h5) file
       whose convergence needs to be checked.
    Outputs:
        If the file passes some basic (not 100% 
        conclusive) convergence checks, it will be 
        returned in h5py format. If it fails, the
        function will raise a IOError, KeyError,
        or ValueError that can be caught and handled
        elsewhere.
    '''

    import h5py
    import numpy as np

    f = h5py.File(file, 'r')
    _ = f['finished'][()]
    shouldBePresent = f['FSABFlow'][()]
    if np.any(np.isnan(shouldBePresent)):
        raise IOError
    if np.all(f['particleFlux_vm_rN'][()] == 0.0): # Indicates result was not stored
        raise IOError

    return f

def combineAndSort(IVvec, DVarr):

    '''
    Inputs:
        IVvec: List or Numpy array of independent
        variable data.
        DVarr: List, or list of lists, or 1D Numpy
        array, or 2D numpy array of dependent 
        variable data. One of the array dimensions
        must match the length of IVvec.
    Outputs:
        A Numpy array in which the independent variable
        lives in the first column, the dependent
        variable series live in subsequent columns, and
        the data is sorted so that the independent
        variable is strictly increasing
    '''

    import numpy as np

    IVvec = np.array(IVvec)
    DVarr = np.array(DVarr)

    if IVvec.ndim != 1:
        raise IOError('<IVvec> must be one-dimensional.')
    
    compatibleLength = IVvec.shape[0]
    DVdim = DVarr.ndim
    DVshape = DVarr.shape
    
    if  DVdim == 1 and DVshape[0] == compatibleLength:
        pass
    elif DVdim == 2 and DVshape[0] == compatibleLength:
        pass
    elif DVdim == 2 and DVshape[1] == compatibleLength:
        DVarr = DVarr.T
    else:
        raise IOError('Neither of the dimensions of <DVarr> were the same as the length of <IVvec>, or <DVarr> was not 2D.')

    combined = np.column_stack((IVvec, DVarr)) # The IV values will be the first column. The data comes in subsequent columns.
    combined = combined[combined[:, 0].argsort()] # This sorts the data so that IVvec values are strictly increasing

    return combined

def relDiff(n1, n2):

    '''
    Inputs:
        n1 and n2: Ints, floats, or NumPy arrays of
                   the same shape.
    Outputs:
        Relative difference of n1 and n2 (or their
        entries) with the same shape as n1 and n2.
    '''

    import numpy as np

    num = n1 - n2
    denom = 0.5 * (n1 + n2)
    
    return np.abs(num / denom)

def thermalVelocity(T, m, units='SI'):

    '''
    Inputs:
        T: Temperature, with units determined by <units>.
        m: Mass, with units determined by <units>.
        units: Defines the unit system used. Options:
               SI: T in joules, m in kilograms.
               keV_mp: T in keV, m in proton masses.
    Ouptuts:
        Thermal velocity in m/s.
    '''

    if units == 'SI':
        pass
    elif units == 'keV_mp':
        from scipy.constants import m_p # Proton mass in kg
        T = T * 1.602176633E-16 # J
        m = m * m_p # kg
    else:
        raise IOError('The requested input unit convention was not recognized.')
    
    return (2 * T / m) ** (1/2)
    
def coulombLog(aDict, bDict):

    '''
    Inputs:
        aDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), and mass (key 'm',
               val in units of proton masses) information
               for one species.
        bDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), and mass (key 'm',
               val in units of proton masses) information
               for another species. If 'z' is -1 in both aDict
               and bDict, the values from aDict alone will be
               used to perform calculations for thermal electron-
               electron collisions. Otherwise, information from
               both dictionaries will be used.
    Outputs:
        An approximation of the Coulomb Logarithm based on
        https://farside.ph.utexas.edu/teaching/plasma/Plasma/node39.html.
        Ion-ion collisions are treated by a formula taken from the
        Neotransp library by Hakan Smith.
    '''

    from scipy.constants import m_p # Proton mass in kg
    from math import log
    
    # Convert units so the formulas are easy to use
    na = aDict['n'] * 1e14 # cm^-3
    nb = bDict['n'] * 1e14 # cm^-3
    ta = aDict['t'] * 1000 # eV
    tb = bDict['t'] * 1000 # eV
    za = aDict['z']
    zb = bDict['z']
    ma_mp = aDict['m']
    mb_mp = bDict['m']

    # Initial check
    if na < 0 or nb < 0 or ta < 0 or tb < 0 or ma_mp < 0 or mb_mp < 0:
        raise IOError('Densities, temperatures, and masses cannot be negative. Please check the inputs.')
    
    # Determine which formula to use and do the calculations
    if za == -1 and zb == -1: # Both species are electrons
        
        if ta < 10:
            return 23 - log(na**(1/2) * ta**(-3/2))
        elif ta > 10:
            return 24 - log(na**(1/2) * ta**(-1))
        else:
            raise ValueError('The algorithm could not determine which formula to use for the Coulomb Logarithm. This sometimes happens in edge cases. Try tweaking the inputs and running again.')
    
    elif za == -1 or zb == -1: # One species is electrons
        
        # Sort out which species is which
        if za == -1:
            te = ta
            ti = tb
            ne = na
            ni = nb
            ze = za
            zi = zb
            me_mp = ma_mp
            mi_mp = mb_mp
        else:
            ti = ta
            te = tb
            ni = na
            ne = nb
            zi = za
            ze = zb
            mi_mp = ma_mp
            me_mp = mb_mp

        # Specify some reused quantities
        me_mi = me_mp / mi_mp
        ti_me_mi = ti * me_mi
        zfac = 10 * zi**2

        # Do the calculations
        if te < ti_me_mi:
            return 30 - log(ne**(1/2) * ti**(-3/2) * zi**(3/2) * mi_mp)
        elif ti_me_mi < te < zfac:
            return 23 - log(ne**(1/2) * zi * te**(-3/2))
        elif ti_me_mi < zfac < te:
            return 24 - log(ne**(1/2) * te**(-1))
        else:
            raise ValueError('The algorithm could not determine which formula to use for the Coulomb Logarithm. This sometimes happens in edge cases. Try tweaking the inputs and running again.')

    else: # Both species are ions

        # Convert some units again for use with a different formula
        ma = ma_mp * m_p # kg
        mb = mb_mp * m_p # kg
        na = na * 1e6 # m^-3
        nb = nb * 1e6 # m^-3

        # Calculations
        t1 = za * zb * (ma + mb)
        t2 = ma * tb + mb * ta
        t3a = na * za**2 / ta
        t3b = nb * zb**2 / tb
        t3 = (t3a + t3b)**(1/2)

        return 30.3 - log(t1 / t2 * t3)

def K_ab(aDict, bDict, K):
    
    '''
    Inputs:
        aDict: Dictionary containing mass (key 'm', val
               float in arbitrary units) and temperature
               (key 't', val float in arbitrary units)
               information for one species.
        bDict: Same as aDict, but for a second species.
               The units used for 'm' and 't' must be
               the same in both dictionaries.
        K: Normalized kinetic energy (mv^2/2)/T of the species
           corresponding to aDict. A value of 1 is often used
           because thermal interactions are frequently of interest.
    Outputs:
        The dimesnionless quantity K^(alpha/beta) in
        Beidler et al, Nuclear Fusion 51 (2011) 076001.
    '''
    
    return bDict['m'] / aDict['m'] * aDict['t'] / bDict['t'] * K

def nu0_ab(aDict, bDict): 

    '''
    Inputs:
        aDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), mass (key 'm',
               val in units of proton masses), and velocity
               (key 'v', val in units of m/s) information
               for one species. Specifically, this is the
               'species of interest' (alpha) from
               Beidler et al, Nuclear Fusion 51 (2011) 076001.
               Note that v must be consistent with the
               normalized kinetic energy K.
        bDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), mass (key 'm',
               val in units of proton masses) information
               for another species. Specifically, this is one
               'background species' (beta) from Beidler et al.
    Outputs:
        The reference collision frequency nu_0^(alpha/beta) in 
        Beidler et al.
    '''

    from scipy.constants import e # Elementary charge in Coulombs
    from scipy.constants import pi
    from scipy.constants import epsilon_0 # In SI units (F/m)
    from scipy.constants import m_p # Proton mass in kg

    # Get inputs in SI units
    nbeta = bDict['n'] * 1e20
    qalpha = aDict['z'] * e
    qbeta = bDict['z'] * e
    coulLog = coulombLog(aDict, bDict)
    malpha = aDict['m'] * m_p
    v = aDict['v']
    
    # Safety checks
    if nbeta <= 0 or malpha <= 0 or v <= 0:
        raise IOError('Species "a" mass, species "a" velocity, and species "b" density must be positive.')

    # Perform calculations
    num = nbeta * (qalpha * qbeta)**2 * coulLog
    denom = 4 * pi * (epsilon_0 * malpha)**2 * v**3

    return num / denom

def nu_ab(aDict, bDict, K):
    
    '''
    Inputs:
        aDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), and mass (key 'm',
               val in units of proton masses) information
               for one species. Specifically, this is the
               'species of interest' (alpha) from
               Beidler et al, Nuclear Fusion 51 (2011) 076001.
        bDict: Dictionary containing density (key 'n', 
               val in units of 10^20 m^-3), temperature
               (key 't', val in units of keV), charge number
               (key 'z', val unitless), and mass (key 'm',
               val in units of proton masses) information
               for another species. Specifically, this is one
               'background species' (beta) from Beidler et al.
        K: Normalized kinetic energy (mv^2/2)/T of the first species.
           A value of 1 is commmon because thermal interactions
           are frequently of interest. Because the implemented
           Coulomb Logarithm calculation is for thermal collisions,
           keeping K=1 is recommended.
    Outputs:
        The collision frequency nu^(alpha/beta) in 
        Beidler et al, Nuclear Fusion 51 (2011) 076001.
    '''

    from scipy.constants import m_e # Electron mass in kg
    from scipy.constants import m_p # Proton mass in kg
    from scipy.constants import pi
    from scipy.special import erf
    from math import exp
    
    # Add velocity key to aDict
    aDict['v'] = thermalVelocity(aDict['t'] * K, aDict['m'], units='keV_mp')

    # Assemble terms
    refColFreq = nu0_ab(aDict, bDict)
    k = K_ab(aDict, bDict, K)
    t1 = erf(k**(1/2)) * (1 - 1/(2 * k))
    t2 = (pi * k)**(-1/2) * exp(-k)

    # Calculate output
    return refColFreq * (t1 + t2)
