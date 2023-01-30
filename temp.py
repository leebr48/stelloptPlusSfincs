# FIXME explain what this script does and its limitations both here and in the eventual args.

# Load necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, sproot, splev, splder
from scipy.integrate import trapezoid

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import combineAndSort
from IO import messagePrinter
from sfincsOutputLib import sfincsRadialAndErScan

#FIXME consider raising warnings for a given subdirectory rather than errors

# Defaults
outdir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/' # FIXME generalize
ErSearchTol = 1.0e-12 # Maximum Jr - this is also used in writeNamelist.py #FIXME might make this an option. Keep in mind you should make it an option for writeNamelist too, and you will perhaps need to pass it to ambipolarSolve in this script so that Sfincs actually performs root finding when you ask it to.

# Locally useful functions
def constructBSpline(dataMat): # Just here to ensure the spline order and smoothing settings are consistent throughout this script
    tck = splrep(dataMat[:,0], dataMat[:,1], k=3, s=0)
    return tck

def findRoots(dataMat, numRootEst, xScan):
    tck = constructBSpline(dataMat)
    estRoots = sproot(tck, mest=numRootEst)
    yEst = splev(xScan, tck)
    return tck, estRoots, yEst

def getErJrData(dataMat, negative=True):
    ErRange = [dataMat[0,0], dataMat[-1,0]]
    extendedErRange = [ErRange[0] - 0.1*np.abs(ErRange[0]), ErRange[1] + 0.1*np.abs(ErRange[1])]

    if negative:
        ErsData = combined[np.where(combined[:,0] < 0)]
        ErScan = np.linspace(extendedErRange[0], 0, num=100)
        numRoots = 1
    else:
        ErsData = combined[np.where(combined[:,0] > 0)]
        ErScan = np.linspace(0, extendedErRange[1], num=100)
        numRoots = 2

    tck, estRoots, JrEst = findRoots(ErsData, numRoots, ErScan)

    return tck, estRoots, ErScan, JrEst

def getAllRootInfo(dataMat, knownRoots):
    negKnownRoots = np.array(knownRoots)[np.where(knownRoots < 0)]
    posKnownRoots = np.array(knownRoots)[np.where(knownRoots > 0)]

    # Negative roots first
    negTck, negEstRoots, negErScan, negJrEst = getErJrData(dataMat, negative=True)
    negStableRoots = determineRootStability(negTck, negKnownRoots)

    # Now positive roots
    posTck, posEstRoots, posErScan, posJrEst = getErJrData(dataMat, negative=False)
    posStableRoots = determineRootStability(posTck, posKnownRoots)
    
    # Now combine everything
    tcks = [negTck, posTck]
    estRoots = np.append(negEstRoots, posEstRoots)
    stableRoots = np.append(negStableRoots, posStableRoots)
    ErScan = np.append(negErScan, posErScan)
    JrScan = np.append(negJrEst, posJrEst)

    return tcks, estRoots, stableRoots, ErScan, JrScan # Note that estRoots will contain any actual roots since there is no smoothing of the polynomials

def findUniqueRoots(actualRoots, estRoots):
    uniqueInds = np.isin(estRoots, actualRoots, invert=True)
    uniqueGuesses = estRoots[uniqueInds]

    return uniqueGuesses

def determineRootStability(tckData, knownRoots):
    if len(knownRoots) != 0:
        tckDer = splder(tckData, n=1)
        ders = splev(knownRoots, tckDer)
        stableInds = np.where(ders > 0) # dJr/dEr > 0 -> stable root
        stableRoots = knownRoots[stableInds]
    else:
        stableRoots = np.array([])
    
    return stableRoots

# Load tools from external library
ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObjCopy', verbose=0) #FIXME generalize address
ErQ = ds.Ersearch(ErQuantity='Er', verbose=0, launch='no') # This has the estimates for the field zeros... note that you only seem to get one per radial directory! #FIXME use proper electric field quantitiy





# Determine where the satisfactory roots are for each radial subdirectory
rootsToUse = [] # FIXME perhaps only save this if you get an answer for each radius? Or create a table so it's clear if there was an error?
for radInd in range(ds.Nradii):
    # Load and sort data from the given radial directory
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    numActualRoots = len(rootErs)
    combined = combineAndSort(ErVals, JrVals)

    # Interpolate between the available data points to determine root stability and guess the position of as-yet-unfound roots
    tcks, estRoots, stableRoots, ErScan, JrScan = getAllRootInfo(combined, rootErs)
    numStableRoots = len(stableRoots)
    uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
    numUniqueRootGuesses = len(uniqueRootGuesses)
    
    # Plot data for interpretation later (if needed)
    plt.figure()
    plt.axhline(y=0, color='black', linestyle='-')
    plt.plot(ErScan, JrScan)
    plt.scatter(combined[:,0], combined[:,1])
    for root in estRoots: #FIXME should you be using estRoots?
        plt.axvline(x=root, color='black', linestyle=':')
    plt.xlabel('$\mathrm{E_r}$')
    plt.ylabel('$\mathrm{J_r}$')
    plt.savefig(outdir+'testImg'+str(radInd)+'.pdf', bbox_inches='tight', dpi=400) #FIXME generalize address

    # Determine if new runs should be launched, or the data processed as-is
    if numActualRoots == 0 and numUniqueRootGuesses == 0:
        errMsg = 'No root could be identified for rN={}. This likely means not enough data was available.'.format(ds.Erscans[radInd].rN)
        errMsg += ' Please generate data for more electric field values and run this script again.'
        messagePrinter(errMsg)
        continue
    elif numActualRoots == 1 and numUniqueRootGuesses == 0: # The identified root is probably the only one
        rootsToUse.append(rootErs[0])
    else: #FIXME this loop stuff... can probably do this at the end of all the checks
        for root in uniqueRootGuesses:
            closestInd = np.argmin(np.abs(root - ErVals))
            ds.Erscans[radInd].launchRun('Er', root, 'nearest', closestInd, ambipolarSolve=True, launchCommand='sbatch') #FIXME generalize 'Er' and 'ambipolarSolve' and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind) #FIXME consider taking away the asking permission to launch, or at least providing an option for it (probably just launch by default)
    # FIXME FIXME FIXME
quit()












# Determine where the satisfactory roots are
rootsToUse = [] # FIXME perhaps only save this if you get an answer for each radius? Or create a table so it's clear if there was an error?
for radInd in range(ds.Nradii):
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    numActualRoots = len(rootErs)
    combined = combineAndSort(ErVals, JrVals)
    
    if numActualRoots in (0,1): #FIXME this needs a lot of testing...
        tcks, estRoots, stableRoots, ErScan, JrScan = getAllRootInfo(combined)
        uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
        numUniqueRootGuesses = len(uniqueRootGuesses)
        if numUniqueRootGuesses == 0 and numActualRoots == 0:
            errMsg = 'No root could be identified for rN={}. This likely means not enough data was available.'.format(ds.Erscans[radInd].rN)
            errMsg += ' Please generate data for more electric field values and run this script again.'
            messagePrinter(errMsg)
            continue
        elif numUniqueRootGuesses == 0 and numActualRoots == 1: # The identified root is probably the only one
            rootsToUse.append(rootErs[0])
        else:
            for root in uniqueRootGuesses:
                closestInd = np.argmin(np.abs(root - ErVals))
                ds.Erscans[radInd].launchRun('Er', root, 'nearest', closestInd, ambipolarSolve=True, launchCommand='sbatch') #FIXME generalize 'Er' and 'ambipolarSolve' and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind) #FIXME consider taking away the asking permission to launch, or at least providing an option for it (probably just launch by default)
    
    elif numActualRoots in (2,3): #FIXME this needs a lot of testing...
        numStableRoots = len(stableRoots)
        if numStableRoots == 1 and numActualRoots == 2:
            tcks, estRoots, stableRoots, ErScan, JrScan = getAllRootInfo(combined)
            uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
            numUniqueRootGuesses = len(uniqueRootGuesses)
            if numUniqueRootGuesses == 0:
                errMsg = 'One stable and one unstable root were found for the rN={} subdirectory, but the location of a second stable root could not be estimated.'.format(ds.Erscans[radInd].rN)
                errMsg += ' Please generate data for more electric field values and run this script again.'
                raise ValueError(errMsg)
            else:
                for root in uniqueRootGuesses:
                    closestInd = np.argmin(np.abs(root - ErVals))
                    ds.Erscans[radInd].launchRun('Er', root, 'nearest', closestInd, ambipolarSolve=True, launchCommand='sbatch') #FIXME generalize 'Er' and 'ambipolarSolve' and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind) #FIXME consider taking away the asking permission to launch, or at least providing an option for it (probably just launch by default)
        elif numStableRoots == 1 and numActualRoots == 3:
            raise ValueError('Three roots were found in the rN={} subdirectory, both two of them appear to be unstable. Something is wrong.'.format(ds.Erscans[radInd].rN))
        elif numStableRoots == 2:
            ionRoot = np.min(stableRoots)
            electronRoot = np.max(stableRoots)
            truncCombined = combined[np.where((ionRoot <= combined[:,0]) & (combined[:,0] <= electronRoot))]
            intVal = trapezoid(truncCombined[:,1], x=truncCombined[:,0]) # FIXME could in principle use polynomial integration, but you'd have to treat the middle bit separately
            if intVal > 0:
                rootsToUse.append(ionRoot) 
            elif intVal < 0:
                rootsToUse.append(electronRoot)
            else:
                raise ValueError('An integral of Jr with respect to Er was identically zero. Something is wrong.')
        else:
            raise ValueError('{} roots were found in the rN={} subdirectory, but {} of them appear to be stable. Something is wrong.'.format(numActualRoots, ds.Erscans[radInd].rN, numStableRoots))

        
    elif numActualRoots > 3:
        raise ValueError('More than three roots were detected, which should not be possible. The search tolerance used to determine whether or not a root exists may need to be smaller.') #FIXME mention an input name if you end up making this an input.
