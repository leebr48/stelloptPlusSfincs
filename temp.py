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
ErSearchTol = 1.0e-8 # Maximum Jr - this is also used in writeNamelist.py #FIXME should be 1.0e-12

# Locally useful function
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
    JrScan = np.append(negJrScan, posJrScan)

    return estRoots, ErScan, JrScan

# Load tools from external library
ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObj', verbose=0) #FIXME generalize address
ErQ = ds.Ersearch(ErQuantity='Er', verbose=0, launch='no') # This has the estimates for the field zeros... note that you only seem to get one per radial directory! #FIXME use proper electric field quantitiy

#FIXME note that if only one root is present, you don't need to do the integral and such - that root is the only answer. But you should still run ambipolarSolve for it!
# Determine where the satisfactory roots are # FIXME all this is busted
for radInd in range(ds.Nradii):
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    numRoots = len(rootErs)
    combined = combineAndSort(ErVals, JrVals)
    if numRoots == 0:
        estRoots, ErScan, JrScan = getAllRootInfo(combined)
        #FIXME fit splines and launch runs at estimated root(s)
        #FIXME using external library launch bits may be easiest...
    elif numRoots == 1:
        #FIXME fit splines and calculate roots. Check if any guesses are unique. If they are, send those to Sfincs with ambipolarSolve. If not, take the one root as being the only one.
    elif numRoots in (2,3):
        #FIXME use Turkin et al to determine which root is correct. Note that you will need to reject the center (unstable) root if it exists. Should also probably check the sign of the roots to be sure you actually have electron and ion roots and not just some artifacts.
    elif numRoots > 3:
        #FIXME this is an artifact of using polynomials. Pull out the three roots nearest Er=0, then use Turkin et al.


quit()

# FIXME that outside library is quite the pain... maybe just fit a polynomial to the data and use it to estimate your zeros? If you do this, might ditch the outside library (if not too arduous).
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
