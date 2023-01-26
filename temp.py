# FIXME explain what this script does and its limitations both here and in the eventual args.

# Load necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, sproot, splev

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import combineAndSort
from sfincsOutputLib import sfincsRadialAndErScan

# Defaults
outdir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/' # FIXME generalize
ErSearchTol = 1.0e-12 # Maximum Jr - this is also used in writeNamelist.py #FIXME might make this an option. Keep in mind you should make it an option for writeNamelist too, and you will perhaps need to pass it to ambipolarSolve in this script so that Sfincs actually performs root finding when you ask it to.

# Locally useful functions
def findRoots(dataMat, numRootEst, xScan):
    tck = splrep(dataMat[:,0], dataMat[:,1], k=3, s=0)
    estRoots = sproot(tck, mest=numRootEst)
    yEst = splev(xScan, tck)
    return estRoots, yEst

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

    EstRoots, JrEst = findRoots(ErsData, numRoots, ErScan)

    return EstRoots, ErScan, JrEst

def getAllRootInfo(dataMat):
    negEstRoots, negErScan, negJrEst = getErJrData(dataMat, negative=True)
    posEstRoots, posErScan, posJrEst = getErJrData(dataMat, negative=False)
    
    estRoots = np.append(negEstRoots, posEstRoots)
    ErScan = np.append(negErScan, posErScan)
    JrScan = np.append(negJrEst, posJrEst)

    return estRoots, ErScan, JrScan

def findUniqueRoots(actualRoots, estRoots):
    uniqueInds = np.isin(estRoots, actualRoots, invert=True)
    uniqueGuesses = estRoots[uniqueInds]

    return uniqueGuesses

# Load tools from external library
ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObj', verbose=0) #FIXME generalize address
ErQ = ds.Ersearch(ErQuantity='Er', verbose=0, launch='no') # This has the estimates for the field zeros... note that you only seem to get one per radial directory! #FIXME use proper electric field quantitiy

# Determine where the satisfactory roots are
for radInd in range(ds.Nradii):
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    numRoots = len(rootErs)
    combined = combineAndSort(ErVals, JrVals)
    if numRoots == 0:
        estRoots, ErScan, JrScan = getAllRootInfo(combined)
        for root in estRoots:
            closestInd = np.argmin(np.abs(root - ErVals))
            ds.Erscans[radInd].launchRun('Er', root, 'nearest', closestInd, ambipolarSolve=True, launchCommand='sbatch') #FIXME generalize 'Er' and 'ambipolarSolve' and 'sbatch' if appropriate (keep ambipolarSolve weirdness in mind)

    elif numRoots == 1:
        pass
        estRoots, ErScan, JrScan = getAllRootInfo(combined)
        uniqueRootGuesses = findUniqueRoots(rootErs, estRoots)
        #FIXME fit splines and calculate roots. Check if any guesses are unique. If they are, send those to Sfincs with ambipolarSolve. If not, take the one root as being the only one.
        #FIXME can perhaps merge this with numRoots==0 case?
    
    elif numRoots in (2,3):
        pass
        # FIXME determine if you have a stable or unstable root! If stable, you're good. If unstable, find the third one.
        #FIXME use Turkin et al to determine which root is correct. Note that you will need to reject the center (unstable) root if it exists. Should also probably check the sign of the roots to be sure you actually have electron and ion roots and not just some artifacts.
    
    elif numRoots > 3:
        raise IOError('More than three roots were detected, which should not be possible. The search tolerance used to determine whether or not a root exists may need to be smaller.') #FIXME mention an input name if you end up making this an input.


quit()

# Plot radial current and estimate zeros
Ers = []
Jrs = []
for radInd in range(ds.Nradii):
    
    # Sort and store Er and Jr values for use later
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    combined = combineAndSort(ErVals, JrVals)
    Ers.append(combined[:,0].tolist())
    Jrs.append(combined[:,1].tolist())

    # Estimate what values of Er give Jr=0 using cubic B-splines
    negEstRoots, negErScan, negJrEst = getErJrData(combined, negative=True)
    posEstRoots, posErScan, posJrEst = getErJrData(combined, negative=False)
    estRoots = np.append(negEstRoots,posEstRoots) #FIXME probably store this or something

    # Plot data for interpretation later
    plt.figure()
    plt.plot(negErScan, negJrEst)
    plt.plot(posErScan, posJrEst)
    plt.scatter(combined[:,0], combined[:,1])
    plt.vlines(estRoots, np.min(combined[:,1]), np.max(combined[:,1]), linestyles='dotted')
    plt.xlabel('$\mathrm{E_r}$')
    plt.ylabel('$\mathrm{J_r}$')
    plt.savefig(outdir+str(radInd)+'.pdf', bbox_inches='tight', dpi=400) #FIXME generalize address
