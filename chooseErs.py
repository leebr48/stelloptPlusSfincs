# This script helps the user identify the root(s) of the ambipolar electric field in a given configuration. Note that each flux surface in a device should have either 1 or 3 roots.
# In the former case, the root is identified with confidence by the script. In the latter case, eq. (A2) of Turkin et al., PoP 18, 022505 (2011) is used to determine the correct root.
# (Two roots are theoretically possible, but the second one would be metastable since it would only occur when the point corresponding to dJr/dEr=0 was on the horizontal axis.)
# Note that this method may produce incorrect results - see, for instance, Maassberg et al., PoP 7, 295-311 (2000).
# A normal workflow would be to use run.py for a given configuration to generate many scans over flux surfaces and radial electric field values (WITHOUT using ambipolarSolve), running
# this script repeatedly until all the roots are identified (which may require looking at the plots produced and choosing radial electric field values manually), running this script
# with the <filter> option to copy only the information for the "correct" electric field values to a new directory, and running ploy.py on that directory to get the "correct" system
# behavior. The plots of Jr vs the radial electric field will typically show a spike when the electric field variable is zero. It is highly suggested that the scanned electric field
# values cover much of this upward spike. Ensuring good scan resolution in this region will help the root finding and integration algorithms. Please note that the spike may be extremely
# thin, especially for inner flux surfaces (since we expect Er=0 on the magnetic axis). The roots will likely be difficult to find in this region. Similarly, if you notice that many
# root guesses are close together and the next guess in the region is strictly less than or greater than all the others in that region, consider manually running a case further from
# these guesses to help the root finding algorithm converge faster. If the algorithm still has convergence troubles, SFINCS may need to be run with higher resolution. Note that the
# final electric field should be relatively smooth, except if the system "switches" between the electron and ion roots - this will look like a step change since the equation from
# Turkin et al. assumes the diffusion coefficient D_E = 0. If the electric field is jagged, the script may have chosen the wrong root. In this case, consider switching the value of the
# electric field on that flux surface (in rootsToUse.txt) to the alternative (found in electronRoots.txt or ionRoots.txt). Keep in mind that you may need to substantially increase the
# run time for SFINCS when includePhi1 is turned on. Again, you should NOT use ambipolarSolve with this script - doing so will create complications that can be unpleasant to deal with.
# It is easier and more reliable to simply run this script repeatedly. As a final note, the internal units used in SFINCS and sfincsScan will almost certainly be different than those
# this script uses when making plots.

# Load necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from collections import Counter
from shutil import copy
from warnings import warn

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import combineAndSort, constructBSpline, relDiff, fixOutputUnits
from IO import getChooseErsArgs, getFileInfo, makeDir, findFiles, messagePrinter, prettyDataLabel, saveTimeStampFile
from sfincsOutputLib import sfincsRadialAndErScan

# Get arguments
args = getChooseErsArgs()
        
# Locally useful functions
def findRoots(dataMat, xScan):
    f = PchipInterpolator(dataMat[:,0], dataMat[:,1], extrapolate=False)
    estRoots = f.roots(extrapolate=False)
    yEst = f(xScan)
    return f, estRoots, yEst

def getErJrData(dataMat, negative=True, numInterpPoints=1000):

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
    ErVals = getattr(sfincsScanInstance, electricFieldVar) # In SFINCS internal units
    conversionFactor = fixOutputUnits(electricFieldVar, 1)
    for root in uniqueRootGuesses: # In physical units
        adjustedUnitsRoot = root / conversionFactor
        closestInd = np.argmin(np.abs(adjustedUnitsRoot - ErVals))
        sfincsScanInstance.launchRun(electricFieldVar, adjustedUnitsRoot, 'nearest', closestInd, sendRunToScheduler=(not args.noRun), launchCommand='sbatch')

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
    # Note that if no electric field scan is present, this function will return nonsense. This behavior is caught another way later, though.
    dataFiles = findFiles('sfincsOutput.h5', sfincsDir, raiseError=True) # Note that sfincsScan breaks if you use a different output file name, so the default is hard-coded in
    subdirsFirst = [address.replace(sfincsDir, '') for address in dataFiles]
    radSubdirTitles = [address.split('/')[1] for address in subdirsFirst]
    elecSubdirTitles = [address.split('/')[2] for address in subdirsFirst]
    radSubdirLabels = [title.split('_')[0] for title in radSubdirTitles]
    elecSubdirLabels = [''.join(i for i in label if not i.isdigit()).replace('-','').replace('.','') for label in elecSubdirTitles]
    countRadLabels = Counter(radSubdirLabels)
    countElecLabels = Counter(elecSubdirLabels)
    mostCommonRadLabel = max(countRadLabels, key=countRadLabels.get)
    mostCommonElecLabel = max(countElecLabels, key=countElecLabels.get)

    return mostCommonRadLabel, mostCommonElecLabel

def filterActualRoots(rootErs, rootJrs, diffTol = 0.01):
    
    # This is not efficient, but for our (small) data sets it should be fine
    rootErsToKeep = np.array([])
    rootJrsToKeep = np.array([])
    for ind1, rootEr1 in enumerate(rootErs):
        with np.errstate(invalid='ignore'): # NumPy will throw a warning before evaluating frac, which makes output confusing
            diffs = relDiff(rootEr1, rootErs) # Compare one root with all the others
        tooSimilarInds = np.where(diffs < diffTol)
        if tooSimilarInds[0].size == 0: # No similar roots -> no comparison needed
            ErsToCompare = rootErs
            JrsToCompare = rootJrs
            bestJrInd = ind1
        else:
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
    if len(ErQuantity) != len(Er): # I/O problem
        return np.nan
    zeroInd = np.where(ErQuantity == 0)
    if not np.all(Er[zeroInd] == 0): # Physics problem
        return np.nan
    newErQuantity = np.delete(ErQuantity, zeroInd)
    newEr = np.delete(Er, zeroInd)
    if len(newErQuantity) != len(newEr): # Physics problem
        return np.nan
    trueCount = np.count_nonzero(newErQuantity * newEr > 0) # All entries True if same sign, all entries False otherwise
    if (trueCount != len(newErQuantity)) and (trueCount != 0): # Physics problem
        return np.nan
    if trueCount != 0:
        return True # ErQuantity and Er have the same sign
    else:
        return False # ErQuantity and Er have different signs

def evaluateIntegral(negF, posF, lowerBound, upperBound):

    if (lowerBound <= 0) and (upperBound <= 0):
        return negF.integrate(lowerBound, upperBound, extrapolate=False)

    elif (lowerBound >= 0) and (upperBound >= 0):
        return posF.integrate(lowerBound, upperBound, extrapolate=False)

    elif (lowerBound <= 0) and (upperBound >= 0):
        return negF.integrate(lowerBound, 0, extrapolate=True) + posF.integrate(0, upperBound, extrapolate=True)
        # Jr is forced to be evaluated very near Er = 0, so extrapolation should be fine here.

    elif (lowerBound >= 0) and (upperBound <= 0):
        return posF.integrate(lowerBound, 0, extrapolate=True) + negF.integrate(0, upperBound, extrapolate=True)
        # Jr is forced to be evaluated very near Er = 0, so extrapolation should be fine here.

    else:
        return np.nan

def findPlotMinMax(inData, margin):
    dataMin = np.min(inData)
    dataMax = np.max(inData)
    dataRange = dataMax - dataMin
    useMin = dataMin - margin * dataRange
    useMax = dataMax + margin * dataRange

    return useMin, useMax

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

# Check how the input directory is organized
radLabel, electricFieldLabel = determineLabels(inDir)

# Load tools from external library
ds = sfincsRadialAndErScan(inDir, verbose=0, ErDefForJr=electricFieldLabel)

# Initial check
if len(ds.Erscans) == 0:
    msg = 'It appears that there are no electric field subdirectories in this SFINCS directory. '
    msg += 'Therefore, this script cannot be used.'
    raise IOError(msg)

# Do work
if not args.filter:
    
    # Determine where the satisfactory roots are for each radial subdirectory
    rootsToUse = [] 
    ionRoots = []
    unstableRoots = []
    electronRoots = []
    integralVals = []
    soloRoots = []
    allRootsLists = [rootsToUse, ionRoots, unstableRoots, electronRoots, integralVals, soloRoots] # Also contains integral values... these are not 'roots', but closely related
    for radInd in range(ds.Nradii):
        
        # Load and sort data from the given radial directory
        dataContainer = ds.Erscans[radInd]
        radVal = getattr(dataContainer, radLabel)[0]
        if dataContainer.includePhi1:
            distr = 'vd'
        else:
            distr = 'vm'
        if electricFieldLabel == 'Er':
            coord = 'rHat'
        else:
            coord = electricFieldLabel.split('d')[-1]
        particleFluxVar = 'particleFlux_'+distr+'_'+coord
        radialCurrentVar = 'radialCurrent_'+distr+'_'+coord
        ErVals = fixOutputUnits(electricFieldLabel, getattr(dataContainer, electricFieldLabel)) # Note this doesn't literally have to be Er, it could be various derivatives of the electric potential
        closeToZero = np.isclose(ErVals, 0, rtol=0, atol=args.zeroErTol[0])
        if True not in closeToZero:
            msg = 'No zero-electric-field case (or anything sufficiently close) was found in the {} = {} subdirectory, '.format(radLabel, radVal)
            msg += 'so this subdirectory will be skipped.'
            messagePrinter(msg)
            recordNoEr(allRootsLists)
            continue
        if 0 in dataContainer.Zs:
            msg = 'A zero-charge particle was found in the output of the {} = {} subdirectory, '.format(radLabel, radVal)
            msg += 'so this subdirectory will be skipped.'
            messagePrinter(msg)
            recordNoEr(allRootsLists)
            continue
        actualEr = fixOutputUnits('Er', getattr(dataContainer, 'Er')) # This is just for comparison purposes
        ErQuantityHasSameSignAsEr = determineErQuantitySign(ErVals, actualEr)
        if np.isnan(ErQuantityHasSameSignAsEr):
            msg = 'For {} = {}, the radial electric field as represented by {} did not seem to be consistent '.format(radLabel, radVal, electricFieldLabel)
            msg += 'with the standard definition. It may have been zero in odd places or not followed the proper sign convention. Please investigate. '
            msg += 'In the meantime, this subdirectory will be skipped.'
            messagePrinter(msg)
            recordNoEr(allRootsLists)
            continue
        JrVals = fixOutputUnits(radialCurrentVar, dataContainer.Jr)
        exactlyZeroCurrentInds = np.where(np.array(JrVals) == 0)
        exactlyZeroErVals = np.array(ErVals)[exactlyZeroCurrentInds]
        if len(exactlyZeroErVals) != 0:
            msg = 'For {} = {}, it appears that the following electric field subdirectories contained a run with zero radial current:\n'.format(radLabel, radVal)
            msg += str(exactlyZeroErVals) + '\n'
            msg += 'This indicates a SFINCS error. Please check and fix these run(s) to get reliable results.\n'
            if not args.allowZeroJr:
                msg += 'In the meantime, this subdirectory will be skipped because it will break the root finding algorithm.'
                messagePrinter(msg)
                recordNoEr(allRootsLists)
                continue
            else:
                msg += 'Because <allowZeroJr> was used, calculations for this subdirectory will continue even though the results will likely be incorrect.'
                messagePrinter(msg)
        unsortedParticleFluxes = fixOutputUnits(particleFluxVar, getattr(dataContainer, particleFluxVar))
        ErParticleFluxes = combineAndSort(ErVals, unsortedParticleFluxes)
        
        # Sort out roots
        rootInds = np.where(np.abs(np.array(JrVals)) <= args.maxRootJr[0])
        allRootErs = np.array(ErVals)[rootInds] # Could contain (effective) duplicates in rare cases
        allRootJrs = np.array(JrVals)[rootInds]
        rootErs, rootJrs = filterActualRoots(allRootErs, allRootJrs)
        numActualRoots = len(rootErs)
        ErJrVals = combineAndSort(ErVals, JrVals)
        if args.print:
            msg = 'For {} = {}, the radial electric field (or proxy) values are:\n'.format(radLabel, radVal)
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

        # Plot and save total current data for interpretation later
        plt.figure()
        plt.axhline(y=0, color='black', linestyle='-', zorder=0)
        plt.scatter(ErJrVals[:,0], ErJrVals[:,1], zorder=5)
        for ErPart, JrPart in zip(ErScan, JrScan):
            plt.plot(ErPart, JrPart, color='tab:blue', zorder=10)
        for root in estRoots:
            if root in rootErs:
                ls = '-'
                if root in stableRoots:
                    clr = 'green'
                    lbl = 'Stable Root'
                else:
                    clr = 'red'
                    lbl = 'Unstable Root'
            else:
                ls = ':'
                clr = 'black'
                lbl = 'Root Guess'
            plt.axvline(x=root, color=clr, linestyle=ls, label=lbl, zorder=15)
        useMin, useMax = findPlotMinMax(ErJrVals[:,1], args.marg[0])
        plt.ylim(bottom=useMin, top=useMax)
        plt.legend(loc='best')
        plt.xlabel(prettyDataLabel(electricFieldLabel))
        plt.ylabel(prettyDataLabel(radialCurrentVar))
        nameBits = basename(inDir) + '-' + radLabel + '_' + str(radVal) + '-' + 'Jr-vs-' + electricFieldLabel
        plotName = nameBits + '.pdf'
        dataName = nameBits + '.dat'
        plt.savefig(join(outDir, plotName), bbox_inches='tight', dpi=400)
        np.savetxt(join(outDir, dataName), ErJrVals)
        plt.close()

        # Also plot individual species flux data
        plt.figure()
        plt.axhline(y=0, color='black', linestyle='-', zorder=0)
        for i in range(1, ErParticleFluxes.shape[1]):
            lbl = 'Spec. ' + str(i)
            plt.plot(ErParticleFluxes[:,0], ErParticleFluxes[:,i], label=lbl, zorder=5)
        useMin, useMax = findPlotMinMax(ErParticleFluxes[:,1:], args.marg[0])
        plt.ylim(bottom=useMin, top=useMax)
        plt.legend(loc='best')
        plt.xlabel(prettyDataLabel(electricFieldLabel))
        plt.ylabel(prettyDataLabel(particleFluxVar))
        nameBits = basename(inDir) + '-' + radLabel + '_' + str(radVal) + '-' + 'particleFluxes-vs-' + electricFieldLabel
        plotName = nameBits + '.pdf'
        dataName = nameBits + '.dat'
        plt.tight_layout()
        plt.savefig(join(outDir, plotName), bbox_inches='tight', dpi=400)
        np.savetxt(join(outDir, dataName), ErParticleFluxes)
        plt.close()
        
        # Also plot grouped (ion and electron) flux data
        plt.figure()
        plt.axhline(y=0, color='black', linestyle='-', zorder=0)
        eInds = np.where(dataContainer.Zs < 0)
        iInds = np.where(dataContainer.Zs > 0)
        if eInds[0].size > 1:
            msg = 'Multiple negative-charge species were detected in the output of the {} = {} subdirectory. '.format(radLabel, radVal)
            msg += 'This is unusual. They will all be labeled "electrons" in the grouped flux plots.'
            warn(msg)
        eFlux = np.sum(ErParticleFluxes[:,1:].T[eInds], axis=0)
        iFlux = np.sum(ErParticleFluxes[:,1:].T[iInds], axis=0)
        plt.plot(ErParticleFluxes[:,0], eFlux, label='Electrons', zorder=5)
        plt.plot(ErParticleFluxes[:,0], iFlux, label='Ions', zorder=6)
        eiFlux = np.column_stack((eFlux, iFlux))
        useMin, useMax = findPlotMinMax(eiFlux, args.marg[0])
        plt.ylim(bottom=useMin, top=useMax)
        plt.legend(loc='best')
        plt.xlabel(prettyDataLabel(electricFieldLabel))
        plt.ylabel(prettyDataLabel(particleFluxVar))
        nameBits = basename(inDir) + '-' + radLabel + '_' + str(radVal) + '-' + 'groupedParticleFluxes-vs-' + electricFieldLabel
        plotName = nameBits + '.pdf'
        dataName = nameBits + '.dat'
        plt.tight_layout()
        plt.savefig(join(outDir, plotName), bbox_inches='tight', dpi=400)
        np.savetxt(join(outDir, dataName), np.column_stack((ErParticleFluxes[:,0], eiFlux)))
        plt.close()

        # Determine if new runs should be launched, or the data processed as-is
        if numActualRoots == 0:
            
            if numUniqueRootGuesses == 0: # No root (real or estimated) can be identified - this is a problem
                printMoreRunsMessage('No root could be identified for {} = {}.'.format(radLabel, radVal))
                recordNoEr(allRootsLists)
            
            else: # A root guess has been identified - launch a run for it
                launchNewRuns(uniqueRootGuesses, dataContainer, electricFieldLabel)
                recordNoEr(allRootsLists)
       
        elif numActualRoots == 1:
            
            if numUniqueRootGuesses == 0: # The fit polynomials could not find any roots beyond the one already identified
                
                if numStableRoots == 0: # The only root that can be found or guessed is unstable - this is a problem
                    printMoreRunsMessage('Only a single, unstable root could be identified for {} = {}.'.format(radLabel, radVal))
                    recordNoEr(allRootsLists)

                else: # The identified root is stable, and no other guesses are apparent - assume the identified root is the only one
                    soloRoots.append(rootErs[0])
                    rootsToUse.append(rootErs[0])
                    ionRoots.append(np.nan)
                    unstableRoots.append(np.nan)
                    electronRoots.append(np.nan)
                    integralVals.append(np.nan)

            else: # If there are other guesses for roots, they should be investigated
                launchNewRuns(uniqueRootGuesses, dataContainer, electricFieldLabel)
                recordNoEr(allRootsLists)

        elif numActualRoots == 2:

            if numStableRoots == 0: # Impossible - something is wrong
                printMoreRunsMessage('Two unstable roots were identified for {} = {}, which should not be possible.'.format(radLabel, radVal))
                recordNoEr(allRootsLists)

            elif numStableRoots == 1: # The other stable root has not been found yet

                if numUniqueRootGuesses == 0: # This is a problem
                    printMoreRunsMessage('Only one stable and one unstable root were identified for {} = {}.'.format(radLabel, radVal))
                    recordNoEr(allRootsLists)

                else: # Investigate the guesses, which will hopefully allow the other stable root to be identified
                    launchNewRuns(uniqueRootGuesses, dataContainer, electricFieldLabel)
                    recordNoEr(allRootsLists)

            else: # Both roots are stable, which likely means the data is incomplete
                
                if numUniqueRootGuesses == 0: # This is a problem
                    printMoreRunsMessage('Two stable roots and no unstable roots were identified for {} = {}.'.format(radLabel, radVal))
                    recordNoEr(allRootsLists)

                else: # Investigate the guesses, which will hopefully allow the unstable root to be identified
                    launchNewRuns(uniqueRootGuesses, dataContainer, electricFieldLabel)
                    recordNoEr(allRootsLists)
                    
        elif numActualRoots == 3:
            
            if numStableRoots in (0, 1, 3): # Does not make sense, something strange is going on
                printMoreRunsMessage('Three roots were identified for {} = {}. {} of these appear to be stable, which should not be possible.'.format(radLabel, radVal, numStableRoots))
                recordNoEr(allRootsLists)

            else: # We have a valid number of total and stable roots
               
                if numUniqueRootGuesses != 0: # Probably a good idea to investigate these guesses just to get better data resolution (since we can only have 1 or 3 roots)
                    launchNewRuns(uniqueRootGuesses, dataContainer, electricFieldLabel)
                    recordNoEr(allRootsLists)

                else: # So far, it looks like we have three roots with the right stability properties
                   
                    # We need to check that the unstable root is between the stable roots
                    unstableRoot = list(set(rootErs) - set(stableRoots))[0]
                    lowerRoot = np.min(stableRoots)
                    upperRoot = np.max(stableRoots)
                    
                    if lowerRoot < unstableRoot < upperRoot: # Everything is in order, so we can choose the correct root using eq. (A2) of Turkin et al., PoP 18, 022505 (2011)
                        
                        negF = fs[0]
                        posF = fs[1]
                        intVal = evaluateIntegral(negF, posF, lowerRoot, upperRoot)
                        # Note that the bounds of the integral in Turkin et al. would switch when Er is defined as a derivative of the potential without the negative sign,
                        # so this ^ form of the integral should be correct.
                        
                        if np.isnan(intVal) or intVal == 0:
                            msg = 'For {} = {}, the integral of Jr with respect to Er was {}. '.format(radLabel, radVal, intVal)
                            msg += 'Something is wrong, so this subdirectory will be skipped.'
                            messagePrinter(msg)
                            recordNoEr(allRootsLists)
                            continue

                        if ErQuantityHasSameSignAsEr:
                            ionRoot = lowerRoot
                            electronRoot = upperRoot
                        else:
                            ionRoot = upperRoot
                            electronRoot = lowerRoot
                            
                        ionRoots.append(ionRoot)
                        unstableRoots.append(unstableRoot)
                        electronRoots.append(electronRoot)
                        integralVals.append(intVal)
                        soloRoots.append(np.nan)
                        
                        if intVal > 0:
                            rootsToUse.append(ionRoot)
                        elif intVal < 0:
                            rootsToUse.append(electronRoot)

                    else:
                        printMoreRunsMessage('For {} = {}, three roots were found and two were stable, but the unstable root was not between the stable ones.'.format(radLabel, radVal))
                        recordNoEr(allRootsLists)

        else: # More than three roots should not be possible
            printMoreRunsMessage('More than three roots were identified for {} = {}, which should not be possible.'.format(radLabel, radVal))
            recordNoEr(allRootsLists)

    # Now perform some checks and save the Er information that was found
    assert ds.Nradii == len(rootsToUse), 'The vector to be written in rootsToUse.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(ionRoots), 'The vector to be written in ionRoots.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(unstableRoots), 'The vector to be written in unstableRoots.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(electronRoots), 'The vector to be written in electronRoots.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(integralVals), 'The vector to be written in integralVals.txt was the wrong length. Something is wrong.'
    assert ds.Nradii == len(soloRoots), 'The vector to be written in soloRoots.txt was the wrong length. Something is wrong.'

    np.savetxt(join(outDir, 'rootsToUse.txt'), rootsToUse)
    np.savetxt(join(outDir, 'ionRoots.txt'), ionRoots)
    np.savetxt(join(outDir, 'unstableRoots.txt'), unstableRoots)
    np.savetxt(join(outDir, 'electronRoots.txt'), electronRoots)
    np.savetxt(join(outDir, 'integralVals.txt'), integralVals)
    np.savetxt(join(outDir, 'soloRoots.txt'), soloRoots)

    # Write a log file
    logStr = 'This directory was last auto-analyzed to determine the correct values of the ambipolar radial electric field on:\n'
    saveTimeStampFile(outDir, 'automatedErDeterminationLog', logStr)
    
    # Closing message
    messagePrinter('The root-choosing algorithms have run on {}.'.format(inDir))
    messagePrinter('Please check the outputs in {} to see the status of the calculations.'.format(outDir))

else:
    
    # Load correct electric field quantities
    loadedErQuantities = np.loadtxt(join(inDir, 'determineEr/rootsToUse.txt'), ndmin=1)
    assert np.any(np.isnan(loadedErQuantities)) == False, 'At least one of the flux surfaces in <sfincsDir> did not have a "correct" electric field specified in rootsToUse.txt.'

    # Copy the correct directories
    for radInd in range(ds.Nradii):
        
        dataContainer = ds.Erscans[radInd]

        ErVals = getattr(dataContainer, electricFieldLabel)
        correctErQuantity = loadedErQuantities[radInd]

        correctInd = np.argmin((ErVals - correctErQuantity) ** 2)
        correctErDir = dataContainer.DataDirs[correctInd]
        
        dirToCopyFrom = join(dataContainer.mainDir, correctErDir)
        dirToCopyTo = dirname(dirToCopyFrom.replace(inDir, outDir)) # Using dirname gets rid of the electric field subdirectory - it is no longer needed
        _ = makeDir(dirToCopyTo)
        copy(join(dirToCopyFrom, 'sfincsOutput.h5'), join(dirToCopyTo, 'sfincsOutput.h5')) # Note that this script has file overwrite powers!
    
    # Write a log file
    logStr = 'This directory was created by copying the SFINCS runs with the "correct" values for the radial electric field from {} on:\n'.format(inDir)
    saveTimeStampFile(outDir, 'automatedCopyLog', logStr)

    # Closing message
    messagePrinter('The data for the "correct" radial electric field has been copied from {} to {}.'.format(inDir, outDir))
