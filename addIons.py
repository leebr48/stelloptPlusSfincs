# This script adds new ions to plasma profiles. The outputs are in STELLOPT format. Input profiles should be quasineutral (and this is checked).
# The electron density profile will not be changed, but already-present ion density profiles may be. Also, the S knot values for the electron
# density profile will be used for all other profiles. Please note that for the purposes of determining the pressure profile, it is assumed that
# the temperature profile for all the ion species matches that of the first ion species.

# Import necessary modules
import numpy as np
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from scipy.linalg import solve

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getAddIonsArgs, getFileInfo, cleanStrings, listifyBEAMS3DFile, makeProfileNames, extractProfileData, extractScalarData, makeStringForStellopt, messagePrinter
from dataProc import nonlinearInterp, relDiff

# Important constant
eVToJ = 1.602176634e-19 # In STELLOPT, temperatures are written in eV but pressures are written in Pa

# Sort out inputs
args = getAddIonsArgs()
inFile, _, _, _, _ = getFileInfo(args.profilesIn[0], '/arbitrary/path', 'arbitrary')
newMasses = args.ms # kg
newCharges = args.zs
newFracs = args.fs

# Handy function
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
profileString = makeStringForStellopt('NI_AUX_S', sVec)
for ionInd in range(totalNumIons):
    fortranInd = ionInd + 1
    profileString += makeStringForStellopt('NI_AUX_F({},:)'.format(fortranInd), dens[:, ionInd])
profileString += makeStringForStellopt(massInName, ionMs)
profileString += makeStringForStellopt(chargeInName, ionZs, integer=True)

presString = makeStringForStellopt('PMASS_TYPE', "'akima_spline'")
presString += makeStringForStellopt('PRES_SCALE', [1])
presString += makeStringForStellopt('AM_AUX_S', sVec)
presString += makeStringForStellopt('AM_AUX_F', pres)

# Print outputs
presMsg = 'The pressure information (relevant for VMEC) is:\n'
presMsg += presString
messagePrinter(presMsg)

profMsg = 'The species information (relevant for BEAMS3D) is:\n'
profMsg += profileString
messagePrinter(profMsg)

messagePrinter('Please replace the appropriate lines in the namelist with the ones printed above.')
