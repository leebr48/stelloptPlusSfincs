# FIXME this script adds new ions to plasma profiles. The outputs are in STELLOPT format.
# FIXME inputs should be quasineutral # FIXME maybe check that?
# FIXME probably note that ne is fixed

# Import necessary modules
import numpy as np
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from scipy.linalg import solve #FIXME needed?

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import cleanStrings, listifyBEAMS3DFile, makeProfileNames, extractProfileData, extractScalarData
from dataProc import nonlinearInterp

# Important constant
eVToJ = 1.602176634e-19 # STELLOPT temperatures are written in eV

# Desired new ion parameters # FIXME use args!
He_frac = 0.05 # Fraction of total ion number density attributed to this species #FIXME kill anything with this variable
He_mass = 6.64647907E-27 # kg #FIXME kill anything with this variable
He_charge = 2 # charge number #FIXME kill anything with this variable
newCharges = [3.0, 4.0] # FIXME pull from args
newFracs = [0.05, 0.10] # FIXME pull from args

# Handy function
def cleanup(ar, integer=False):
    if integer:
        precision = 0
        ar = ar.astype(int)
    else:
        precision = 10
    newstr = np.array2string(ar, separator='     ', precision=precision, max_line_width=100).replace('[','').replace(']','')
    return newstr

# Retrieve data
inFile = '/cobra/u/lebra/data/w7x/reactor/input.W7X_REACTOR_woptim_forSfincs' # FIXME use args
varsToFind = ['NE', 'NI', 'TE', 'TI']
prefixesOfInterest = cleanStrings(varsToFind)
listifiedInFile = listifyBEAMS3DFile(inFile)
profileVarsOfInterest = makeProfileNames(prefixesOfInterest)
profileData = extractProfileData(listifiedInFile, profileVarsOfInterest)
ders = {}
for key,val in profileData.items():
    ders[key] = 0
profileFitFuncs = nonlinearInterp(profileData, ders, pchip=True) # Use in case the S vectors are not uniform
scalarVarsOfInterest = cleanStrings(['NI_AUX_M', 'NI_AUX_Z']) # STELLOPT has the mass and charge of electrons built in, so only the ions need to be specified
scalarData = extractScalarData(listifiedInFile, scalarVarsOfInterest)

# Calculations
sVec = profileData['ne']['iv'][0] # FIXME generalize, or just note that this is what's done?
numOldIons = len(profileFitFuncs['ni'])
numNewIons = len(newCharges)
ionZs = scalarData['z'] + newCharges # this sets the species order
totalNumIons = len(ionZs)
allZs = [-1] + ionZs
pres = []
for sVal in sVec: # Calculations must be done one flux surface at a time
    
    # Declare terms of the equation Ax=b
    A = []
    b = []
    
    # Charge conservation
    A.append(ionZs)
    localNE = float(profileFitFuncs['ne'][0](sVal))
    b.append(localNE)

    # New ion particle balance
    for newIonInd, newIonFrac in enumerate(newFracs): # one equation for each new ion
        coeffs = [newIonFrac] * totalNumIons
        coeffs[numOldIons + newIonInd] = newIonFrac - 1
        A.append(coeffs)
        b.append(0)

    if numOldIons > 1: # FIXME be sure to test with and without!
        # Particle balance closure - redundant, but necessary to get a square matrix
        for oldIonInd, oldIonFunc in enumerate(profileFitFuncs['ni'][1:]):
            r = oldIonFunc(sVal) / profileFitFuncs['ni'][0](sVal)
            coeffs = [0] * totalNumIons
            coeffs[0] = r
            coeffs[oldIonInd + 1] = -1
            A.append(coeffs)
            b.append(0)

    # Solve the system for the densities
    x = solve(A, b) # FIXME note that the accuracy of the solution (like when you do np.dot(A, x) == b) isn't amazing, but maybe it's close enough?

    # Calculate the pressure on the given flux surface
    ne_te = localNE * profileFitFuncs['te'][0](sVal)
    ni_ti = x * profileFitFuncs['ti'][0](sVal) # FIXME assumes that the temperature of all ions is the same
    pres.append((ne_te + np.sum(ni_ti)) * eVToJ)

quit()
# Calculations # FIXME caluclate pressure! # FIXME multi-species?
nD = ne / 2 / (1 + 2 * He_frac / (1 - He_frac))
nT = nD
nHe = 2 * He_frac * nD / (1 - He_frac)

mi = np.append(mi, He_mass)
zi = np.append(zi, He_charge)

# Output #FIXME multi-species?
print('  NI_AUX_F(1,:) =', cleanup(nD))
print('  NI_AUX_F(2,:) =', cleanup(nT))
print('  NI_AUX_F(3,:) =', cleanup(nHe))
print('  NI_AUX_M =', cleanup(mi))
print('  NI_AUX_Z =', cleanup(zi, integer=True))

print('******************') # FIXME should these checks be asserts instead?
QN = -1 * ne + zi[0] * nD + zi[1] * nT + zi[2] * nHe
HeIonDensFrac = nHe / (nD + nT + nHe)
print('Quasineutrality check: these should (in principle) be zero, but you must keep numerical error in mind: ', QN)
print('Accuracy check: these should be the impurity fraction: ', HeIonDensFrac)
# FIXME maybe print to a file by default and have a print option in args?
# FIXME don't forget to write/print AM_AUX_S with PMASS_TYPE = 'akima_spline' and PRES_SCALE =  1.00000000000000E+00
