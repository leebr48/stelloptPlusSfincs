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

# Load tools from external library
ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObj', verbose=0) #FIXME generalize address
ErQ = ds.Ersearch(ErQuantity='Er', verbose=0, launch='no') # This has the estimates for the field zeros... note that you only seem to get one per radial directory! #FIXME use proper electric field quantitiy

# FIXME that outside library is quite the pain... maybe just fit a polynomial to the data and use it to estimate your zeros? If you do this, might ditch the outside library (if not too arduous).
# FIXME can use scipy.interpolate.splrep and scipy.interpolate.sproot
# Plot radial current and estimate zeros
Ers = []
Jrs = []
for radInd in range(ds.Nradii):
    ErVals = ds.Erscans[radInd].Er
    JrVals = ds.Erscans[radInd].Jr
    combined = combineAndSort(ErVals, JrVals)
    Ers.append(combined[:,0].tolist())
    Jrs.append(combined[:,1].tolist())
    ErRange = [combined[0,0], combined[-1,0]]
    tck = splrep(combined[:,0], combined[:,1], k=3) #FIXME polynomials are NOT working... maybe exponentials? Excel thinks those don't work either... What does neotransp use?
    estRoots = sproot(tck, mest=3) #FIXME must have at least 8 knots to use this
    print(estRoots)
    ErScan = np.linspace(ErRange[0] - 0.1*np.abs(ErRange[0]), ErRange[1] + 0.1*np.abs(ErRange[1]))
    JrEst = splev(ErScan, tck)
    plt.figure()
    plt.plot(ErScan, JrEst)
    plt.scatter(combined[:,0], combined[:,1])
    plt.xlabel('$\mathrm{E_r}$')
    plt.ylabel('$\mathrm{J_r}$')
    plt.savefig(outdir+str(radInd)+'.pdf', bbox_inches='tight', dpi=400)

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
