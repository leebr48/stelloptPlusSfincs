# FIXME Given pressure and toroidal flux information from a VMEC wout_* file and bootstrap information from a SFINCS directory,
# this script will calculate the "correct" bootstrap current for the given magnetic configuration. This could be used, for
# instance, to refine an optimized configuration until the bootstrap is self-consistent. The information is output in STELLOPT
# format so it can be used with VMEC.
# This script numerically integrates equation (17) in "Computing vmec’s ac current profile and curtor from a bootstrap current code"
# by Matt Landreman, which is included in the STELLOPT directory in SHARE/doc/computing_vmec_AC_profile_from_a_bootstrap_current_code.

#FIXME note that this script will not work with an Er scan
#FIXME note that extrapolation is used!
#FIXME note that you need sfincs runs throughtout the volume to get a decent result!

# FIXME perhaps test on self-consistent bootstrap configuration?

# Load necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
from scipy.interpolate import PchipInterpolator
from scipy.io.netcdf import netcdf_file
from scipy.constants import mu_0
from scipy.integrate import odeint
from numpy import linspace, pi

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from sfincsOutputLib import sfincsRadialAndErScan, sfincsScan
from dataProc import fixOutputUnits

# Sort out inputs # FIXME command line!
sfincsDir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/seventhObj2_d_correctEr'
woutFile = '/u/lebra/src/stelloptPlusSfincs/outsideTest/seventhObj2_d/wout_W7X_REACTOR_woptim_forSfincs.00126.nc'

# Iniital check
test = sfincsRadialAndErScan(sfincsDir, verbose=0)

if len(test.Erscans) != 0:
    msg = 'It appears that there are electric field subdirectories in this SFINCS directory. '
    msg += 'Therefore, this script cannot be used.'
    raise IOError(msg)

# Load SFINCS information
ds = sfincsScan(sfincsDir, verbose=0)

psiN = ds.psiN
FSABjHat = fixOutputUnits('FSABjHat', ds.FSABjHat)
FSABHat2 = fixOutputUnits('FSABHat2', ds.FSABHat2)

f_jB = PchipInterpolator(psiN, FSABjHat) #NOTE: extrapolation is allowed since SFINCS is usually not run for rN=0 or rN=1
f_B2 = PchipInterpolator(psiN, FSABHat2)

# Load VMEC information
woutFile = netcdf_file(woutFile, mmap=False)

signgs = woutFile.variables['signgs'][()]
ns = woutFile.variables['ns'][()]
pres = woutFile.variables['pres'][()]
dpsids = woutFile.variables['phips'][()][-1] # Any valid index except 0 will return the same (correct) value

vmec_s = linspace(0, 1, num=ns)
f_p = PchipInterpolator(vmec_s, pres)
f_dpds = f_p.derivative()

# Integrate
f_dIds = lambda I, s: (2 * pi * dpsids * f_jB(s) - mu_0 * I * f_dpds(s)) / f_B2(s)
I0 = 0 # no enclosed current on the axis
eval_s = linspace(0, 1, num=10*ns)
Isolved = odeint(f_dIds, I0, eval_s)

# Get output quantities
f_I = PchipInterpolator(eval_s, Isolved)
f_dIds_out = f_I.derivative()
s_for_file = ds.psiN.tolist()
if 0 not in s_for_file:
    s_for_file.insert(0,0)
if 1 not in s_for_file:
    s_for_file.append(1)

curtor = signgs * f_I(1) # equation (24) in Matt's document
outI = f_dIds_out(s_for_file)
# FIXME finish! Can use that makeString function
