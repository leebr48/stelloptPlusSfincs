# FIXME this script adds new ions to plasma profiles. The outputs are in STELLOPT format.
# FIXME inputs should be quasineutral (and this is checked)
# FIXME probably note that ne is fixed
# FIXME note that the S values for NE will be the points on which the pressure profile is evaluated 
# FIXME for the purposes of determining the pressure, the temperature profile for all ion species is assumed to match that of the first ion species. (T same for all species in a reactor anyway)

# FIXME ensure the args have checks on the length of inputs and so on

# Import necessary modules
import numpy as np
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from scipy.linalg import solve

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import cleanStrings, listifyBEAMS3DFile, makeProfileNames, extractProfileData, extractScalarData, messagePrinter
from dataProc import nonlinearInterp, relDiff

# Important constant
eVToJ = 1.602176634e-19 # In STELLOPT, temperatures are written in eV but pressures are written in Pa

# Sort out inputs
inFile = '/raven/u/lebra/src/stelloptPlusSfincs/temp/input.DT' # FIXME use args
newMasses = [6.64647907E-27] # kg # FIXME pull from args
newCharges = [2.0] # FIXME pull from args # FIXME ensure these are ints
newFracs = [0.05] # FIXME pull from args # FIXME ensure the sum of these is less than or equal to 1

# Handy functions
def cleanup(inAr, integer=False):
    ar = np.array(inAr)
    if integer:
        precision = 0
        ar = ar.astype(int)
    else:
        precision = 10 # Standard precision for STELLOPT
    newstr = np.array2string(ar, separator='     ', precision=precision, max_line_width=100).replace('[','').replace(']','')
    return newstr

def makeString(varName, data, integer=False, namePrefix=' '*2, nameSuffix=' = ', eol='\n'):
    if type(data) != str:
        dataStr = cleanup(data, integer=integer)
    else:
        dataStr = data
    return namePrefix + varName + nameSuffix + dataStr + eol

def checkQN(NE, ionZs, NIs, thresh=1e-6, failmsg='Quasineutrality check failed.'):
    QN = -1 * NE + np.dot(ionZs, NIs)
    relQN = np.abs(QN / (NE + np.sum(NIs)))
    assert relQN < thresh, failmsg # thresh is not zero because of numerical error

# Retrieve data
profileVarsToFind = ['NE', 'NI', 'TE', 'TI']
massInName = 'NI_AUX_M'
chargeInName = 'NI_AUX_Z'
scalarVarsToFind = [massInName, chargeInName]
prefixesOfInterest = cleanStrings(profileVarsToFind)
listifiedInFile = listifyBEAMS3DFile(inFile)
profileVarsOfInterest = makeProfileNames(prefixesOfInterest)
profileData = extractProfileData(listifiedInFile, profileVarsOfInterest)
ders = {}
for key,val in profileData.items():
    ders[key] = 0
profileFitFuncs = nonlinearInterp(profileData, ders, pchip=True) # Use in case the S vectors are not uniform
scalarVarsOfInterest = cleanStrings(scalarVarsToFind) # STELLOPT has the mass and charge of electrons built in, so only the ions need to be specified
scalarData = extractScalarData(listifiedInFile, scalarVarsOfInterest)

# Some easy administrative things
ionZs = scalarData['z'] + newCharges # This sets the species order
ionMs = scalarData['m'] + newMasses
numOldIons = len(profileFitFuncs['ni'])
numNewIons = len(newCharges)
totalNumIons = len(ionZs)

# Calculations
sVec = profileData['ne']['iv'][0]
pres = []
for sInd, sVal in enumerate(sVec): # Calculations must be done one flux surface at a time

    # Initial quasineutrality check
    localNE = float(profileFitFuncs['ne'][0](sVal))
    localNIs = [oldIonFunc(sVal) for oldIonFunc in profileFitFuncs['ni']]
    checkQN(localNE, scalarData['z'], localNIs, failmsg='The inputs do not appear to fulfill quasineutrality.')
    
    # Declare terms of the equation Ax=b
    A = []
    b = []
    
    # Charge conservation
    A.append(ionZs)
    b.append(localNE)

    # New ion particle balance
    for newIonInd, newIonFrac in enumerate(newFracs): # one equation for each new ion
        coeffs = [newIonFrac] * totalNumIons
        coeffs[numOldIons + newIonInd] = newIonFrac - 1
        A.append(coeffs)
        b.append(0)

    if numOldIons > 1:
        # Particle balance closure - redundant, but necessary to get a square matrix
        for oldIonInd, oldIonFunc in enumerate(profileFitFuncs['ni'][1:]):
            r = oldIonFunc(sVal) / profileFitFuncs['ni'][0](sVal)
            coeffs = [0] * totalNumIons
            coeffs[0] = r
            coeffs[oldIonInd + 1] = -1
            A.append(coeffs)
            b.append(0)

    # Solve the system for the densities
    x = solve(A, b)

    # Record the densities
    if sInd == 0:
        dens = x
    else:
        dens = np.vstack((dens, x))

    # Calculate the pressure on the given flux surface
    ne_te = localNE * profileFitFuncs['te'][0](sVal)
    ni_ti = x * profileFitFuncs['ti'][0](sVal)
    pres.append((ne_te + np.sum(ni_ti)) * eVToJ)

    # Check results
    checkQN(localNE, ionZs, x, failmsg='The outputs do not appear to fulfill quasineutrality. The matrix solve likely went wrong.')
    for ionInd in range(numOldIons, totalNumIons):
        nIonUnderConsideration = x[ionInd]
        achievedFrac = nIonUnderConsideration / np.sum(x)
        desiredFrac = newFracs[ionInd - numOldIons]
        diff = relDiff(achievedFrac, desiredFrac)
        assert diff <= 1e-6, 'The achieved fraction(s) for the added ion(s) did not match the desired value(s). The matrix solve likely went wrong.' # diff is not zero because of numerical error

# Perform final checks, just in case
assert dens.shape[0] == len(pres), 'The density and pressure arrays have inconsistent dimensions. Something is wrong.'
assert dens.shape[1] == totalNumIons, 'The ion density array does not seem to contain the correct number of species. Something is wrong.'

# Prepare outputs
profileString = makeString('NI_AUX_S', sVec)
for ionInd in range(totalNumIons):
    fortranInd = ionInd + 1
    profileString += makeString('NI_AUX_F({},:)'.format(fortranInd), dens[:, ionInd])
profileString += makeString(massInName, ionMs)
profileString += makeString(chargeInName, ionZs, integer=True)

presString = makeString('PMASS_TYPE', "'akima_spline'")
presString += makeString('PRES_SCALE', [1])
presString += makeString('AM_AUX_S', sVec)
presString += makeString('AM_AUX_F', pres)

# Print outputs
presMsg = 'The pressure information (relevant for VMEC) is:\n'
presMsg += presString
messagePrinter(presMsg)

profMsg = 'The profile and species information (relevant for BEAMS3D) is:\n'
profMsg += profileString
messagePrinter(profMsg)

messagePrinter('Note that you MUST delete the old versions of these lines from the namelist before adding the new ones!')
