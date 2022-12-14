from booz_xform import Booz_xform

woutAddress = '/raven/u/lebra/src/stelloptPlusSfincs/RHSmode/2x2Part2/wout_W7X_EIM_n15_e30_i16_8SH2.nc' # FIXME make general, of course
surface_inds = [2, 15, 40, 61] # FIXME should be set by input - you need to coordinate with other codes
booz_toroidal_harmonics = 51 # FIXME how much resolution do you need? - for test case, 101 was much slower and seemed to give identical results
booz_poloidal_harmonics = 51 # FIXME how much resolution do you need? - for test case, 101 was much slower and seemed to give identical results

# Initialize the class
b = Booz_xform()

# Load our equilibrium and set the resolution parameters
b.read_wout(woutAddress) # Note we could also use .read_boozmn(boozmn_*.nc) if desired
b.nboz = booz_toroidal_harmonics
b.mboz = booz_poloidal_harmonics

# Set the surfaces (by their indices) on which we'd like to do calculations
b.compute_surfs = surface_inds

# Do the calculations
b.run()

# Extract the stuff we need
Gs = b.Boozer_G # Only on compute surfaces
Is = b.Boozer_I # Only on compute surfaces
iotas = [b.iota[ind] for ind in b.compute_surfs] # Only on compute surfaces

outDic = {}
for surface_ind, G, I, iota in zip(surface_inds, Gs, Is, iotas):
    outDic[surface_ind] = {'G':G, 'I':I, 'iota':iota}

from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from os import environ
from subprocess import run

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))

from dataProc import getBoozerInformation

outDic = getBoozerInformation(woutAddress, surface_inds)
print(outDic)
