# Given pressure and toroidal flux information from a VMEC wout_* file and bootstrap information from a SFINCS directory,
# this script will calculate the "correct" bootstrap current for the given magnetic configuration. This could be used, for
# instance, to refine an optimized configuration until the bootstrap is self-consistent. The information is output in STELLOPT
# format so it can be used with VMEC. Please note that to get an accurate picture of the bootstrap current inside a stellarator,
# one must formally integrate neoclassically-derived quantities from the magnetic axis to the last closed flux surface. SFINCS
# should therefore be run over as much of the volume as possible. Polynomial extrapolation will be used to "fill in the ends" 
# since SFINCS is not typically run too close to the magnetic axis, and is often not run too close to the last closed flux surface
# either.
# This script numerically integrates equation (17) in "Computing vmec’s ac current profile and curtor from a bootstrap current code"
# by Matt Landreman, which is included in the STELLOPT directory in SHARE/doc/computing_vmec_AC_profile_from_a_bootstrap_current_code.
# The results of this script should be fairly accurate - they are within a few percent of the results in https://doi.org/10.1063/5.0098166,
# for instance, but the source of the discrepancy is not known.

# Load necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys
from scipy.interpolate import PchipInterpolator
from scipy.io import netcdf_file
from scipy.constants import mu_0
from scipy.integrate import odeint
from numpy import linspace, pi

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from sfincsOutputLib import sfincsRadialAndErScan, sfincsScan
from dataProc import fixOutputUnits, createVMECGrids
from IO import getBootstrapArgs, getFileInfo, makeStringForStellopt, messagePrinter

# Sort out inputs
args = getBootstrapArgs()
woutFile, _, _, _, _ = getFileInfo(args.eqIn[0], '/arbitrary/path', 'arbitrary')
_, _, _, sfincsDir, _ = getFileInfo('arbitrary', args.sfincsDir[0], 'arbitrary')

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
pres = woutFile.variables['presf'][()] # Pressure on the full grid - should be the "real" value (half grid is used for derivatives)
dpsids = woutFile.variables['phips'][()][-1] # Any valid index except 0 will return the same (correct) value

_, fullgrid = createVMECGrids(ns)
f_p = PchipInterpolator(fullgrid, pres)
f_dpds = f_p.derivative()

# Integrate
f_dIds = lambda I, s: (2 * pi * dpsids * f_jB(s) - mu_0 * I * f_dpds(s)) / f_B2(s)
I0 = 0 # no enclosed current on the axis
eval_s = linspace(0, 1, num=100*ns)
Isolved = odeint(f_dIds, I0, eval_s)

# Get output quantities
f_I = PchipInterpolator(eval_s, signgs * Isolved)
s_for_file = ds.psiN.tolist()
if 0 not in s_for_file:
    s_for_file.insert(0, 0)
if 1 not in s_for_file:
    s_for_file.append(1)

curtor = f_I(1) # equation (24) in Matt's document, but signgs is already built in
outI = f_I(s_for_file).flatten()

# Print output quantities
curString = 'The current information (relevant for VMEC) is:\n'
curString += makeStringForStellopt('NCURR', [1])
curString += makeStringForStellopt('CURTOR', curtor)
curString += makeStringForStellopt('PCURR_TYPE', "'akima_spline_I'")
curString += makeStringForStellopt('AC_AUX_S', s_for_file)
curString += makeStringForStellopt('AC_AUX_F', outI)
messagePrinter(curString)
messagePrinter('Please replace the appropriate lines in the namelist with the ones printed above.')
