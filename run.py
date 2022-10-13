# This script is the top-level wrapper from which all the other scripts can be called.
# You can trigger scripts to be written and run from here.

# Import necessary modules
from subprocess import run
from IO import getArgs
import writeProfiles
import writeNamelist
import writeBatch

# Get command line arguments
args = getArgs()

# Write requested files
if not args.noProfiles:
    writeProfiles.run()

if not args.noNamelist:
    writeNamelist.run()

if not args.noBatch:
    writeBatch.run()

# Call sfincsScan
run(['sfincsScan']) 
