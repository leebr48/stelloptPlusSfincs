# This script is the top-level wrapper from which all the other scripts can be called.
# You can trigger scripts to be written and run from here.

# Import necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from os import environ
from subprocess import run

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getRunArgs, adjustInputLengths, makeDir, messagePrinter
import writeProfiles
import writeNamelist
import writeBatch

# Get command line arguments
args = getRunArgs()

# Organize the directories that we will work in
inLists = {'profilesIn':args.profilesIn, 'eqIn':args.eqIn, 'saveLoc':args.saveLoc}
IOlists, longLists, maxLen = adjustInputLengths(inLists)

# If saveLoc was not specified, decide whether to use the profiles or equilibria locations
if all([item == None for item in IOlists['saveLoc']]):
    if 'profilesIn' in longLists:
        saveDefaultTarget = IOlists['profilesIn']
    else:
        saveDefaultTarget = IOlists['eqIn']
else:
    saveDefaultTarget = [join(location,'arbitrary') for location in IOlists['saveLoc']]

# Loop through the working directories
for i in range(maxLen):
    profilesInUse = IOlists['profilesIn'][i]
    eqInUse = IOlists['eqIn'][i]
    actualSaveLoc = dirname(saveDefaultTarget[i])

    # Make target directory if it does not exist
    outDir = makeDir(actualSaveLoc) # Note that this script has file overwrite powers!

    # Write requested files
    if not args.noProfiles:
        writeProfiles.run(profilesInUse, outDir)

    if not args.noNamelist:
        writeNamelist.run(profilesInUse, outDir, eqInUse)

    if not args.noBatch:
        writeBatch.run(profilesInUse, outDir)

    # Call sfincsScan if requested
    if not args.noRun:
        execLoc = join(environ['SFINCS_PATH'],'fortran/version3/utils/sfincsScan')

        if args.noConfirm:
            run([execLoc, 'arbitraryCommandLineArg'], cwd=outDir)
        else:
            run([execLoc], cwd=outDir)

    messagePrinter('Setup and/or run task(s) {} of {} completed in {}.'.format(i+1, maxLen, outDir))
