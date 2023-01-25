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
ErSearchTol = 1.0e-12 # Maximum Jr - this is also used in writeNamelist.py

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

# Load tools from external library
ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObj', verbose=0) #FIXME generalize address
ErQ = ds.Ersearch(ErQuantity='Er', verbose=0, launch='no') # This has the estimates for the field zeros... note that you only seem to get one per radial directory! #FIXME use proper electric field quantitiy

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

quit()

#FIXME note that if only one root is present, you don't need to do the integral and such - that root is the only answer. But you should still run ambipolarSolve for it!
# Determine where the satisfactory roots are # FIXME all this is busted
ionRootErs = []
unstableRootErs = []
electronRootErs = []
for radInd in range(ds.Nradii):
    ErVals = Ers[radInd]
    JrVals = Jrs[radInd]
    rootInds = np.where(np.abs(np.array(JrVals)) <= ErSearchTol)
    rootErs = np.array(ErVals)[rootInds]
    if len(rootErs) == 3:
       ionRootErs.append(rootErs[0])
       unstableRootErs.append(rootErs[1])
       electronRootErs.append(rootErs[2])
    else:
       print('Only {} root(s) have been identified for rN = {}.'.format(len(rootErs), ds.Erscans[radInd].rN[0])) 
