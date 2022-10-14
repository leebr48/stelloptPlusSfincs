# This script is the top-level wrapper from which all the other scripts can be called.
# You can trigger scripts to be written and run from here.

# Import necessary modules
from subprocess import run
from IO import getArgs, getFileInfo
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

# Call sfincsScan if requested
if not args.noRun:
    _, _, _, outFilePath, _ = getFileInfo(args.profilesIn[0], args.saveLoc[0], 'arbitrary')
    
    if args.noConfirm:
        run(['sfincsScan', 'arbitraryCommandLineArg'], cwd=outFilePath)
    else:
        run(['sfincsScan'], cwd=outFilePath) 
