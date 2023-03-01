# FIXME Given pressure and toroidal flux information from a VMEC wout_* file and bootstrap information from a SFINCS directory,
# this script will calculate the "correct" bootstrap current for the given magnetic configuration. This could be used, for
# instance, to refine an optimized configuration until the bootstrap is self-consistent. The information is output in STELLOPT
# format so it can be used with VMEC.

# Load necessary modules
from os.path import dirname, abspath, join, basename
from inspect import getfile, currentframe
import sys

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from sfincsOutputLib import sfincsScan

# Sort out inputs # FIXME command line!
inDir = '/u/lebra/src/stelloptPlusSfincs/outsideTest/seventhObj2_d_correctEr'

# Load SFINCS information
ds = sfincsScan(inDir, verbose=0) # FIXME the script won't work with an Er scan (probably good)... but make sure you know how to use this
print(dir(ds))
print(ds.rN)
