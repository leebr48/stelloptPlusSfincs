# FIXME explain what this script does and its limitations both here and in the eventual args. You should also explain that some manual checking must be done with the output plots to ensure all the roots have been found.
#FIXME note that Er should be pretty smooth... not sure if you want to check that somehow in the script, or just tell users to check it by eye

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
from IO import makeDir, messagePrinter
from sfincsOutputLib import sfincsRadialAndErScan

#FIXME try to generalize the radial labels
#FIXME what will you actually run once you have the correct Er's? Your phi1 script, but modified??

# Defaults
indir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObjCopy' # FIXME generalize
outdir = makeDir(join(indir, 'determineEr')) # FIXME generalize... if necessary
ErSearchTol = 1.0e-12 # Maximum Jr - this is also used in writeNamelist.py #FIXME might make this an option. Keep in mind you should make it an option for writeNamelist too, and you will perhaps need to pass it to ambipolarSolve in this script so that Sfincs actually performs root finding when you ask it to.

# Locally useful bits
def constructBSpline(dataMat): # Just here to ensure the spline order and smoothing settings are consistent throughout this script
    tck = splrep(dataMat[:,0], dataMat[:,1], k=3, s=0)
    return tck

def findRoots(dataMat, numRootEst, xScan):
    tck = constructBSpline(dataMat)
    estRoots = sproot(tck, mest=numRootEst) # Should not extrapolate outside of provided data
    yEst = splev(xScan, tck)
    return tck, estRoots, yEst

def getErJrData(dataMat, negative=True):
    ErRange = [dataMat[0,0], dataMat[-1,0]]
    extendedErRange = [ErRange[0] - 0.1*np.abs(ErRange[0]), ErRange[1] + 0.1*np.abs(ErRange[1])]

    if negative:
        ErsData = ErJrVals[np.where(ErJrVals[:,0] < 0)]
        ErScan = np.linspace(extendedErRange[0], 0, num=100)
        numRoots = 1
    else:
        ErsData = ErJrVals[np.where(ErJrVals[:,0] > 0)]
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

def launchNewRuns(uniqueRootGuesses, sfincsScanInstance):
    ErVals = sfincsScanInstance.Er
    for root in uniqueRootGuesses:
        closestInd = np.argmin(np.abs(root - ErVals))
        sfincsScanInstance.launchRun('Er', root, 'nearest', closestInd, ambipolarSolve=True, launchCommand='sbatch') #FIXME generalize 'Er' and 'ambipolarSolve' and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind) #FIXME consider taking away the asking permission to launch, or at least providing an option for it (probably just launch by default)
        # FIXME depending on how the permissions work, you may need to print a notification when you auto-launch a run

def printMoreRunsMessage(customString):
    standardLittleDataErrorMsg = ' This likely means not enough data was available.'
    standardLittleDataErrorMsg += ' Please check the plot generated for this flux surface,'
    standardLittleDataErrorMsg += ' generate data for more electric field values, and run this script again.'

    customString += standardLittleDataErrorMsg
    messagePrinter(customString)

def recordNoEr(listOfLists):
    for singleList in listOfLists:
        singleList.append(np.nan)

# Load tools from external library
ds = sfincsRadialAndErScan(indir, verbose=0)

# Determine where the satisfactory roots are for each radial subdirectory
rootsToUse = [] 
ionRoots = []
electronRoots = []
soloRoots = []
allRootsLists = [rootsToUse, ionRoots, electronRoots, soloRoots]
for radInd in range(ds.Nradii):
    # Load and sort data from the given radial directory
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    numActualRoots = len(rootErs)
    ErJrVals = combineAndSort(ErVals, JrVals)

    # Interpolate between the available data points to determine root stability and guess the position of as-yet-unfound roots
    tcks, estRoots, stableRoots, ErScan, JrScan = getAllRootInfo(ErJrVals, rootErs)
    numStableRoots = len(stableRoots)
    uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
    numUniqueRootGuesses = len(uniqueRootGuesses)
    
    # Plot data for interpretation later (if needed)
    plt.figure()
    plt.axhline(y=0, color='black', linestyle='-')
    plt.plot(ErScan, JrScan)
    plt.scatter(ErJrVals[:,0], ErJrVals[:,1])
    for root in estRoots:
        plt.axvline(x=root, color='black', linestyle=':')
    for root in rootErs:
        plt.axvline(x=root, color='black', linestyle='-')
    plt.xlabel('$\mathrm{E_r}$') # FIXME generalize label - I think you have a function for that?
    plt.ylabel('$\mathrm{J_r}$')
    plt.savefig(join(outdir, 'testImg'+str(radInd)+'.pdf'), bbox_inches='tight', dpi=400) #FIXME generalize address

    # Determine if new runs should be launched, or the data processed as-is
    if numActualRoots == 0:
        
        if numUniqueRootGuesses == 0: # No root (real or estimated) can be identified - this is a problem
            printMoreRunsMessage('No root could be identified for rN={}.'.format(ds.Erscans[radInd].rN))
            recordNoEr(allRootsLists)
            continue
        
        else: # A root guess has been identified - launch a run for it
            launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd])
            recordNoEr(allRootsLists)
   
    elif numActualRoots == 1:
        
        if numUniqueRootGuesses == 0: # The fit polynomials could not find any roots beyond the one already identified
            
            if numStableRoots == 0: # The only root that can be found or guessed is unstable - this is a problem
                printMoreRunsMessage('Only a single, unstable root could be identified for rN={}.'.format(ds.Erscans[radInd].rN))
                recordNoEr(allRootsLists)
                continue

            else: # The identified root is stable, and no other guesses are apparent - assume the identified root is the only one
                soloRoots.append(rootErs[0])
                rootsToUse.append(rootErs[0])

        else: # If there are other guesses for roots, they should be investigated
            launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd])
            recordNoEr(allRootsLists)

    elif numActualRoots == 2:

        if numStableRoots == 0: # Impossible - something is wrong
            printMoreRunsMessage('Two unstable roots were identified for rN={}, which should not be possible.'.format(ds.Erscans[radInd].rN))
            recordNoEr(allRootsLists)
            continue

        elif numStableRoots == 1: # The other stable root has not been found yet

            if numUniqueRootGuesses == 0: # This is a problem
                printMoreRunsMessage('Only one stable and one unstable root were identified for rN={}.'.format(ds.Erscans[radInd].rN))
                recordNoEr(allRootsLists)
                continue

            else: # Investigate the guesses, which will hopefully allow the other stable root to be identified
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd])
                recordNoEr(allRootsLists)

        else: # Both roots are stable, which likely means the data is incomplete
            
            if numUniqueRootGuesses == 0: # This is a problem
                printMoreRunsMessage('Two stable roots and no unstable roots were identified for rN={}.'.format(ds.Erscans[radInd].rN))
                recordNoEr(allRootsLists)
                continue

            else: # Investigate the guesses, which will hopefully allow the unstable root to be identified
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd])
                recordNoEr(allRootsLists)
                
    elif numActualRoots == 3:
        
        if numStableRoots in (0, 1, 3): # Does not make sense, something strange is going on
            printMoreRunsMessage('Three roots were identified for rN={}. {} of these appear to be stable, which should not be possible.'.format(ds.Erscans[radInd].rN, numStableRoots))
            recordNoEr(allRootsLists)
            continue

        else: # We have a valid number of total and stable roots
           
            if numUniqueRootGuesses != 0: # Probably a good idea to investigate these guesses just to get better data resolution (since we can only have 1 or 3 roots)
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd])
                recordNoEr(allRootsLists)

            else: # So far, it looks like we have three roots with the right stability properties
               
                # We need to check that the unstable root is between the stable roots
                unstableRoot = list(set(rootErs) - set(stableRoots))[0]
                ionRoot = np.min(stableRoots)
                electronRoot = np.max(stableRoots)

                if ionRoot < unstableRoot and unstableRoot < electronRoot: # Everything is in order, so we can choose the correct root using eq. (A2) of Turkin et al., PoP 18, 022505 (2011)
                    
                    truncatedErJrVals = ErJrVals[np.where((ionRoot <= ErJrVals[:,0]) & (ErJrVals[:,0] <= electronRoot))]
                    intVal = trapezoid(truncatedErJrVals[:,1], x=truncatedErJrVals[:,0]) # FIXME could in principle use polynomial integration, but you'd have to treat the middle bit separately
                    
                    if intVal == 0:
                        messagePrinter('The integral of Jr with respect to Er was identically zero for rN={}. Something is wrong.'.format(ds.Erscans[radInd].rN))
                        recordNoEr(allRootsLists)
                        continue

                    ionRoots.append(ionRoot)
                    electronRoots.append(electronRoot)
                    
                    if intVal > 0:
                        rootsToUse.append(ionRoot)
                    elif intVal < 0:
                        rootsToUse.append(electronRoot)

                else:
                    printMoreRunsMessage('For rN={}, three roots were found and two were stable, but the unstable root was not between the stable ones.'.format(ds.Erscans[radInd].rN))
                    recordNoEr(allRootsLists)
                    continue

    else: # More than three roots should not be possible
        printMoreRunsMessage('More than three roots were identified for rN={}, which should not be possible.'.format(ds.Erscans[radInd].rN))
        recordNoEr(allRootsLists)
        continue

# Now save the Er information that was found
np.savetxt(join(outdir, 'rootsToUse.txt'), rootsToUse)
np.savetxt(join(outdir, 'ionRoots.txt'), ionRoots)
np.savetxt(join(outdir, 'electronRoots.txt'), electronRoots)
np.savetxt(join(outdir, 'soloRoots.txt'), soloRoots)
