# This script helps the user identify the root(s) of the ambipolar electric field in a given configuration. Note that each flux surface in a device should have either 1 or 3 roots.
# In the former case, the root is identified with confidence by the script. In the latter case, eq. (A2) of Turkin et al., PoP 18, 022505 (2011) is used to determine the correct root.
# A normal workflow would be to use run.py for a given configuration to generate many scans over flux surfaces and radial electric field values (WITHOUT using ambipolarSolve), running
# this script repeatedly until all the roots are identified (which may sometimes require looking at the plots produced and choosing radial electric field values manually), running
# this script with the <filter> option to copy only the information for the "correct" electric field values to a new directory, and running ploy.py on that directory to get the
# "correct" system behavior. The plots of Jr vs the radial electric field will typically show a spike when the electric field variable is zero. It is suggested that the scanned
# electric field values span into part of this upward spike, but do not get too close to the peak. This is because polynomial integration is used for most of the points on this
# plot, but the polynomials cannot fit the full spike properly. So, placing some points near the peak, but not so near that the fit polynomials become inaccurate, should give the
# best results. Note that the final electric field should be relatively smooth, except if the system "switches" between the electron and ion roots - this will look like a step
# change since the equation from Turkin et al. assumes the diffusion coefficient D_E = 0. If the electric field is jagged, the script may have chosen the wrong root. In this case,
# consider switching the value of the electric field on that flux surface (in rootsToUse.txt) to the alternative (found in electronRoots.txt or ionRoots.txt). Keep in mind
# that you may need to substantially modify the run time for SFINCS when includePhi1 is turned on. Also keep in mind that you should NOT use ambipolarSolve with this script - 
# doing so will create complications that can be unpleasant to deal with. It is easier and more reliable to simply run this script repeatedly.

# Load necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from collections import Counter
from shutil import copy

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import combineAndSort, constructBSpline
from IO import getChooseErArgs, getFileInfo, makeDir, findFiles, messagePrinter, prettyDataLabel, saveTimeStampFile
from sfincsOutputLib import sfincsRadialAndErScan

# FIXME lots of testing is needed!
# FIXME the method of Turkin probably requires Er... just force everything to use that?
# FIXME probably mention what the vertical lines on the plots mean?
# FIXME can you plot things after finding the roots and such, so that you don't need to run the code once at the end? (Might already essentially be doing this)
# FIXME note that odd numbers of jobs work well?
# FIXME ensure that run.py <time> isn't zero or negative?

# Get arguments
args = getChooseErArgs()

# Locally useful functions
def findRoots(dataMat, xScan):
    f = PchipInterpolator(dataMat[:,0], dataMat[:,1], extrapolate=False)
    estRoots = f.roots(extrapolate=False)
    yEst = f(xScan)
    return f, estRoots, yEst

def getErJrData(dataMat, negative=True, numInterpPoints=100):

    if negative:
        ErsData = dataMat[np.where(dataMat[:,0] <= 0)]
    else:
        ErsData = dataMat[np.where(dataMat[:,0] >= 0)]
        
    ErScan = np.linspace(np.min(ErsData[:,0]), np.max(ErsData[:,0]), num=numInterpPoints)

    f, estRoots, JrEst = findRoots(ErsData, ErScan)

    return f, estRoots, ErScan, JrEst

def getAllRootInfo(dataMat, knownRoots, ErQuantityHasSameSignAsEr):
    negKnownRoots = np.array(knownRoots)[np.where(knownRoots < 0)]
    posKnownRoots = np.array(knownRoots)[np.where(knownRoots > 0)]

    # Negative roots first
    negF, negEstRoots, negErScan, negJrEst = getErJrData(dataMat, negative=True)
    negStableRoots = determineRootStability(negF, negKnownRoots, ErQuantityHasSameSignAsEr)

    # Now positive roots
    posF, posEstRoots, posErScan, posJrEst = getErJrData(dataMat, negative=False)
    posStableRoots = determineRootStability(posF, posKnownRoots, ErQuantityHasSameSignAsEr)
    
    # Now combine everything
    fs = [negF, posF]
    estRoots = np.append(negEstRoots, posEstRoots)
    stableRoots = np.append(negStableRoots, posStableRoots)
    ErScan = [negErScan, posErScan]
    JrScan = [negJrEst, posJrEst]

    return fs, estRoots, stableRoots, ErScan, JrScan # Note that estRoots will contain any actual roots since there is no smoothing of the polynomials

def findUniqueRoots(actualRoots, estRoots):
    uniqueInds = np.isin(estRoots, actualRoots, invert=True)
    uniqueGuesses = estRoots[uniqueInds]

    return uniqueGuesses

def determineRootStability(f, knownRoots, ErQuantityHasSameSignAsEr):
    if len(knownRoots) != 0:
        fDer = f.derivative()
        ders = fDer(knownRoots)
        if ErQuantityHasSameSignAsEr:
            stableInds = np.where(ders > 0) # dJr/dEr > 0 -> stable root
        else:
            stableInds = np.where(ders < 0) # Using Er definition/representation with opposite sign -> root stability condition flips
        stableRoots = knownRoots[stableInds]
    else:
        stableRoots = np.array([])
    
    return stableRoots

def launchNewRuns(uniqueRootGuesses, sfincsScanInstance, electricFieldVar):
    ErVals = getattr(sfincsScanInstance, electricFieldVar)
    for root in uniqueRootGuesses:
        closestInd = np.argmin(np.abs(root - ErVals))
        sfincsScanInstance.launchRun(electricFieldVar, root, 'nearest', closestInd, sendRunToScheduler=(not args.noRun), launchCommand='sbatch')

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
    # This will only work if the radial directories are named 'var_val', rather than just given an integer number.
    # Directories created with stelloptPlusSfincs will always follow this naming convention.
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

def relDiff(n1, n2):
    num = n1 - n2
    denom = 0.5 * (n1 + n2)
    
    return np.abs(num / denom)

def filterActualRoots(rootErs, rootJrs, diffTol = 0.01): #FIXME might need to change (raise?) diffTol if this doesn't work as expected
    
    # This is not efficient, but for our (small) data sets it should be fine
    rootErsToKeep = np.array([])
    rootJrsToKeep = np.array([])
    for ind1, rootEr1 in enumerate(rootErs):
       diffs = relDiff(rootEr1, rootErs) # Compare one root with all the others
       tooSimilarInds = np.where(diffs < diffTol)
       ErsToCompare = rootErs[tooSimilarInds]
       JrsToCompare = rootJrs[tooSimilarInds]
       bestJrInd = np.argmin(np.abs(JrsToCompare))
       rootErsToKeep = np.append(rootErsToKeep, ErsToCompare[bestJrInd])
       rootJrsToKeep = np.append(rootJrsToKeep, JrsToCompare[bestJrInd])

    sortMat = np.unique(combineAndSort(rootErsToKeep, rootJrsToKeep), axis=0)

    ErsOut = sortMat[:,0]
    JrsOut = sortMat[:,1]

    return ErsOut, JrsOut

def determineErQuantitySign(ErQuantity, Er):
    return np.all(ErQuantity * Er >= 0) # True if same sign, False if different signs

# Sort out directories
_, _, _, inDir, _ = getFileInfo('/arbitrary/path', args.sfincsDir[0], 'arbitrary')

if args.saveLoc[0] is None:
    if not args.filter:
        outDir = join(inDir, 'determineEr')
    else:
        outDir = inDir + '_correctEr'
else:
    _, _, _, outDir, _ = getFileInfo('/arbitrary/path', args.saveLoc[0], 'arbitrary')

_ = makeDir(outDir)

# Load tools from external library
ds = sfincsRadialAndErScan(inDir, verbose=0)

# Check how the input directory is organized
radLabel, electricFieldLabel = determineLabels(inDir)

if not args.filter:
    
    # Determine where the satisfactory roots are for each radial subdirectory
    rootsToUse = [] 
    ionRoots = []
    electronRoots = []
    soloRoots = []
    allRootsLists = [rootsToUse, ionRoots, electronRoots, soloRoots]
    for radInd in range(ds.Nradii):
        
        # Load and sort data from the given radial directory
        ErVals = getattr(ds.Erscans[radInd], electricFieldLabel) # Note this doesn't literally have to be Er, it could be various derivatives of the potential
        actualEr = getattr(ds.Erscans[radInd], 'Er') # This is only for sign comparisons
        ErQuantityHasSameSignAsEr = determineErQuantitySign(ErVals, actualEr)
        JrVals = ds.Erscans[radInd].Jr
        rootInds = np.where(np.abs(np.array(JrVals)) <= args.maxRootJr[0])
        allRootErs = np.array(ErVals)[rootInds] # Could contain (effective) duplicates in rare cases
        allRootJrs = np.array(JrVals)[rootInds]
        rootErs, rootJrs = filterActualRoots(allRootErs, allRootJrs)
        numActualRoots = len(rootErs)
        ErJrVals = combineAndSort(ErVals, JrVals)
        if args.print:
            msg = 'For {} = {}, the radial electric field (or proxy) values are:\n'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0])
            msg += str(ErJrVals[:,0])+'\n'
            msg += 'The corresponding radial current values are:\n'
            msg += str(ErJrVals[:,1])
            messagePrinter(msg)

        # Interpolate between the available data points to determine root stability and guess the position of as-yet-unfound roots
        fs, allEstRoots, stableRoots, ErScan, JrScan = getAllRootInfo(ErJrVals, rootErs, ErQuantityHasSameSignAsEr)
        numStableRoots = len(stableRoots)
        estRoots, _ = filterActualRoots(np.append(rootErs, allEstRoots), np.append(rootJrs, [10 ** 9]*len(allEstRoots))) # Eliminate estimated roots if they are too close to real ones
        uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
        numUniqueRootGuesses = len(uniqueRootGuesses)

        # Plot data for interpretation
        plt.figure()
        plt.axhline(y=0, color='black', linestyle='-')
        plt.scatter(ErJrVals[:,0], ErJrVals[:,1])
        for ErPart, JrPart in zip(ErScan, JrScan):
            plt.plot(ErPart, JrPart, color='tab:blue')
        for root in estRoots:
            if root in rootErs:
                ls = '-'
                if root in stableRoots:
                    clr = 'green'
                else:
                    clr = 'red'
            else:
                ls = ':'
                clr = 'black'
            plt.axvline(x=root, color=clr, linestyle=ls)
        dataMin = np.min(ErJrVals[:,1])
        dataMax = np.max(ErJrVals[:,1])
        dataRange = dataMax - dataMin
        marg = 0.02 # The Matplotlib margins function would not work properly for some reason, so it has to be done manually
        useMin = dataMin - marg * dataRange
        useMax = dataMax + marg * dataRange
        plt.ylim(bottom=useMin, top=useMax)
        plt.xlabel(prettyDataLabel(electricFieldLabel))
        plt.ylabel(prettyDataLabel('radialCurrent_vm_rN')) # vm or vd (no Phi1 or Phi1) shouldn't matter in this case
        plotName = basename(inDir) + '-' + radLabel + '_' + str(getattr(ds.Erscans[radInd], radLabel)[0]) + '-' + 'Jr-vs-' + electricFieldLabel + '.pdf'
        plt.savefig(join(outDir, plotName), bbox_inches='tight', dpi=400)

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
                    lowerRoot = np.min(stableRoots)
                    upperRoot = np.max(stableRoots)

                    if lowerRoot < unstableRoot < upperRoot: # Everything is in order, so we can choose the correct root using eq. (A2) of Turkin et al., PoP 18, 022505 (2011)
                        
                        negF = fs[0]
                        posF = fs[1]

                        # FIXME test this integral stuff a lot, it's important!
                        negIntVal = negF.integrate(np.min(ErJrVals[:,0]), 0, extrapolate=False) # FIXME be sure 0 is included in the domain - makes sense in general since there is a spike there
                        posIntVal = posF.integrate(0, np.max(ErJrVals[:,0]), extrapolate=False)

                        intVal = negIntVal + posIntVal
                        
                        if np.isnan(intVal) or intVal == 0:
                            messagePrinter('The integral of Jr with respect to Er was NaN or identically zero for {} = {}. Something is wrong.'.format(radLabel, getattr(ds.Erscans[radInd], radLabel)[0]))
                            recordNoEr(allRootsLists)
                            continue

                        if ErQuantityHasSameSignAsEr:
                            ionRoot = lowerRoot
                            electronRoot = upperRoot
                        else:
                            ionRoot = upperRoot
                            electronRoot = lowerRoot
                        
                        ionRoots.append(ionRoot)
                        electronRoots.append(electronRoot)
                        soloRoots.append(np.nan)
                        
                        if intVal > 0: # FIXME does this work any more with multiple Er defs?
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

    # Now perform some checks and save the Er information that was found
    assert ds.Nradii == len(rootsToUse), 'The vector to be written in rootsToUse.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(ionRoots), 'The vector to be written in ionRoots.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(electronRoots), 'The vector to be written in electronRoots.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(soloRoots), 'The vector to be written in soloRoots.txt was the wrong length. Something is wrong.'

    np.savetxt(join(outDir, 'rootsToUse.txt'), rootsToUse)
    np.savetxt(join(outDir, 'ionRoots.txt'), ionRoots)
    np.savetxt(join(outDir, 'electronRoots.txt'), electronRoots)
    np.savetxt(join(outDir, 'soloRoots.txt'), soloRoots)

    # Write a log file
    logStr = 'This directory was last auto-analyzed to determine the correct values of the ambipolar radial electric field on:\n'
    saveTimeStampFile(outDir, 'automatedErDeterminationLog', logStr)
    
    # Closing message
    messagePrinter('The root-choosing algorithms have run on {}.'.format(inDir))
    messagePrinter('Please check the outputs in {} to see the status of the calculations.'.format(outDir))

else:
    
    # Load correct electric field quantities
    loadedErQuantities = np.loadtxt(join(inDir, 'determineEr/rootsToUse.txt'))
    assert any(np.isnan(loadedErQuantities)) == False, 'At least one of the flux surfaces in <sfincsDir> did not have a "correct" electric field specified in rootsToUse.txt.'

    # Copy the correct directories
    for radInd in range(ds.Nradii):
        
        ErVals = getattr(ds.Erscans[radInd], electricFieldLabel)
        correctErQuantity = loadedErQuantities[radInd]

        correctInd = np.argmin((ErVals - correctErQuantity) ** 2)
        correctErDir = ds.Erscans[radInd].DataDirs[correctInd]
        
        dirToCopyFrom = join(ds.Erscans[radInd].mainDir, correctErDir)
        dirToCopyTo = dirname(dirToCopyFrom.replace(inDir, outDir)) # Using dirname gets rid of the electric field subdirectory - it is no longer needed
        _ = makeDir(dirToCopyTo)
        copy(join(dirToCopyFrom, 'sfincsOutput.h5'), join(dirToCopyTo, 'sfincsOutput.h5')) # Note that this script has file overwrite powers!
    
    # Write a log file
    logStr = 'This directory was created by copying the SFINCS runs with the "correct" values for the radial electric field from {} on:\n'.format(inDir)
    saveTimeStampFile(outDir, 'automatedCopyLog', logStr)

    # Closing message
    messagePrinter('The data for the "correct" radial electric field has been copied from {} to {}.'.format(inDir, outDir))
