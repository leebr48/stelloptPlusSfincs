# FIXME explain what this script does and its limitations both here and in the eventual args. You should also explain that some manual checking must be done with the output plots to ensure all the roots have been found.
#FIXME note that Er should be pretty smooth (except for electron-ion root change, since you assume D_E=0)... not sure if you want to check that somehow in the script, or just tell users to check it by eye

# Load necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, sproot, splev, splder, splint
from scipy.integrate import trapezoid
from collections import Counter

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import combineAndSort
from IO import makeDir, findFiles, messagePrinter, prettyDataLabel
from sfincsOutputLib import sfincsRadialAndErScan

#FIXME what will you actually run once you have the correct Er's? Your phi1 script, but modified?? (Shouldn't be terrible). What about just using launchRun?
#FIXME lots of testing is needed!

# Administrative matters
indir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObjCopy' # FIXME generalize AND regularize
outdir = makeDir(join(indir, 'determineEr')) # FIXME generalize... if necessary
ErSearchTol = 1.0e-12 # Maximum Jr - this is also used in writeNamelist.py #FIXME might make this an option. Keep in mind you should make it an option for writeNamelist too, and you will perhaps need to pass it to ambipolarSolve in this script so that Sfincs actually performs root finding when you ask it to.

# Locally useful bits
def constructBSpline(dataMat): # Just here to ensure the spline order and smoothing settings are consistent throughout this script #FIXME you might ought to move this elsewhere? Other scripts could perhaps use it...
    tck = splrep(dataMat[:,0], dataMat[:,1], k=3, s=0)
    return tck

def findRoots(dataMat, numRootEst, xScan):
    tck = constructBSpline(dataMat)
    estRoots = sproot(tck, mest=numRootEst) # Should not extrapolate outside of provided data
    yEst = splev(xScan, tck)
    return tck, estRoots, yEst

def getErJrData(dataMat, negative=True, numInterpPoints=100):
    ErRange = [dataMat[0,0], dataMat[-1,0]]

    if negative:
        ErsData = ErJrVals[np.where(ErJrVals[:,0] < 0)]
        ErScanRange = [ErRange[0], 0]
        numRoots = 1
    else:
        ErsData = ErJrVals[np.where(ErJrVals[:,0] > 0)]
        ErScanRange = [0, ErRange[1]]
        numRoots = 2
        
    ErScan = np.linspace(*ErScanRange, num=numInterpPoints)

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

def launchNewRuns(uniqueRootGuesses, sfincsScanInstance, electricFieldVar):
    ErVals = getattr(sfincsScanInstance, electricFieldVar)
    for root in uniqueRootGuesses:
        closestInd = np.argmin(np.abs(root - ErVals))
        sfincsScanInstance.launchRun(electricFieldVar, root, 'nearest', closestInd, ambipolarSolve=True, sendRunToScheduler=True, launchCommand='sbatch') #FIXME generalize 'ambipolarSolve' and the schedulerRun and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind) #FIXME consider taking away the asking permission to launch, or at least providing an option for it (probably just launch by default)
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

def determineLabels(sfincsDir):
    # A lot of this is probably overkill, but it should help add a level of automation safety
    dataFiles = findFiles('sfincsOutput.h5', sfincsDir, raiseError=True) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in
    subdirsFirst = [address.replace(sfincsDir, '') for address in dataFiles]
    radSubdirTitles = [address.split('/')[1] for address in subdirsFirst]
    elecSubdirTitles = [address.split('/')[2] for address in subdirsFirst]
    radSubdirLabels = [title.split('_')[0] for title in radSubdirTitles]
    elecSubdirLabels = [''.join(i for i in label if not i.isdigit()).replace('-','') for label in elecSubdirTitles]
    countRadLabels = Counter(radSubdirLabels)
    countElecLabels = Counter(elecSubdirLabels)
    mostCommonRadLabel = max(countRadLabels, key=countRadLabels.get)
    mostCommonElecLabel = max(countElecLabels, key=countElecLabels.get)

    return mostCommonRadLabel, mostCommonElecLabel

def getSplineBounds(tck):
    return [np.min(tck[0]), np.max(tck[0])]

# Load tools from external library
ds = sfincsRadialAndErScan(indir, verbose=0)

# Check how the input directory is organized
radLabel, electricFieldLabel = determineLabels(indir)

# Determine where the satisfactory roots are for each radial subdirectory
rootsToUse = [] 
ionRoots = []
electronRoots = []
soloRoots = []
allRootsLists = [rootsToUse, ionRoots, electronRoots, soloRoots]
for radInd in range(ds.Nradii):
    
    # Load and sort data from the given radial directory
    ErVals = getattr(ds.Erscans[radInd], electricFieldLabel)
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
    plt.xlabel(prettyDataLabel(electricFieldLabel))
    plt.ylabel(prettyDataLabel('radialCurrent_vm_rN')) # vm or vd (no Phi1 or Phi1) shouldn't matter in this case
    plt.savefig(join(outdir, 'testImg'+str(radInd)+'.pdf'), bbox_inches='tight', dpi=400) #FIXME generalize address

    # Determine if new runs should be launched, or the data processed as-is
    if numActualRoots == 0:
        
        if numUniqueRootGuesses == 0: # No root (real or estimated) can be identified - this is a problem
            printMoreRunsMessage('No root could be identified for {} = {}.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
            recordNoEr(allRootsLists)
            continue
        
        else: # A root guess has been identified - launch a run for it
            launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd], electricFieldLabel)
            recordNoEr(allRootsLists)
   
    elif numActualRoots == 1:
        
        if numUniqueRootGuesses == 0: # The fit polynomials could not find any roots beyond the one already identified
            
            if numStableRoots == 0: # The only root that can be found or guessed is unstable - this is a problem
                printMoreRunsMessage('Only a single, unstable root could be identified for {} = {}.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                recordNoEr(allRootsLists)
                continue

            else: # The identified root is stable, and no other guesses are apparent - assume the identified root is the only one
                soloRoots.append(rootErs[0])
                rootsToUse.append(rootErs[0])
                ionRoots.append(np.nan)
                electronRoots.append(np.nan)

        else: # If there are other guesses for roots, they should be investigated
            launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd], electricFieldLabel)
            recordNoEr(allRootsLists)

    elif numActualRoots == 2:

        if numStableRoots == 0: # Impossible - something is wrong
            printMoreRunsMessage('Two unstable roots were identified for {} = {}, which should not be possible.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
            recordNoEr(allRootsLists)
            continue

        elif numStableRoots == 1: # The other stable root has not been found yet

            if numUniqueRootGuesses == 0: # This is a problem
                printMoreRunsMessage('Only one stable and one unstable root were identified for {} = {}.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                recordNoEr(allRootsLists)
                continue

            else: # Investigate the guesses, which will hopefully allow the other stable root to be identified
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd], electricFieldLabel)
                recordNoEr(allRootsLists)

        else: # Both roots are stable, which likely means the data is incomplete
            
            if numUniqueRootGuesses == 0: # This is a problem
                printMoreRunsMessage('Two stable roots and no unstable roots were identified for {} = {}.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                recordNoEr(allRootsLists)
                continue

            else: # Investigate the guesses, which will hopefully allow the unstable root to be identified
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd], electricFieldLabel)
                recordNoEr(allRootsLists)
                
    elif numActualRoots == 3:
        
        if numStableRoots in (0, 1, 3): # Does not make sense, something strange is going on
            printMoreRunsMessage('Three roots were identified for {} = {}. {} of these appear to be stable, which should not be possible.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0], numStableRoots))
            recordNoEr(allRootsLists)
            continue

        else: # We have a valid number of total and stable roots
           
            if numUniqueRootGuesses != 0: # Probably a good idea to investigate these guesses just to get better data resolution (since we can only have 1 or 3 roots)
                launchNewRuns(uniqueRootGuesses, ds.Erscans[radInd], electricFieldLabel)
                recordNoEr(allRootsLists)

            else: # So far, it looks like we have three roots with the right stability properties
               
                # We need to check that the unstable root is between the stable roots
                unstableRoot = list(set(rootErs) - set(stableRoots))[0]
                ionRoot = np.min(stableRoots)
                electronRoot = np.max(stableRoots)

                if ionRoot < unstableRoot and unstableRoot < electronRoot: # Everything is in order, so we can choose the correct root using eq. (A2) of Turkin et al., PoP 18, 022505 (2011)
                    
                    negTck = tcks[0]
                    posTck = tcks[1]

                    negBounds = getSplineBounds(negTck)
                    posBounds = getSplineBounds(posTck)
                    truncatedErJrVals = ErJrVals[np.where((negBounds[1] <= ErJrVals[:,0]) & (ErJrVals[:,0] <= posBounds[0]))] # Use poly. int. over the fit domain, trapezoidal int. elsewhere

                    negIntVal = splint(*negBounds, negTck)
                    posIntVal = splint(*posBounds, posTck)
                    middleIntVal = trapezoid(truncatedErJrVals[:,1], x=truncatedErJrVals[:,0])

                    intVal = negIntVal + posIntVal + middleIntVal
                    
                    if intVal == 0:
                        messagePrinter('The integral of Jr with respect to Er was identically zero for {} = {}. Something is wrong.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                        recordNoEr(allRootsLists)
                        continue

                    ionRoots.append(ionRoot)
                    electronRoots.append(electronRoot)
                    
                    if intVal > 0:
                        rootsToUse.append(ionRoot)
                    elif intVal < 0:
                        rootsToUse.append(electronRoot)

                else:
                    printMoreRunsMessage('For {} = {}, three roots were found and two were stable, but the unstable root was not between the stable ones.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                    recordNoEr(allRootsLists)
                    continue

    else: # More than three roots should not be possible
        printMoreRunsMessage('More than three roots were identified for {} = {}, which should not be possible.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
        recordNoEr(allRootsLists)
        continue

# Now save the Er information that was found
np.savetxt(join(outdir, 'rootsToUse.txt'), rootsToUse)
np.savetxt(join(outdir, 'ionRoots.txt'), ionRoots)
np.savetxt(join(outdir, 'electronRoots.txt'), electronRoots)
np.savetxt(join(outdir, 'soloRoots.txt'), soloRoots)
