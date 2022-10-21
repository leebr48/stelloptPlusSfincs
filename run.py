# This script is the top-level wrapper from which all the other scripts can be called.
# You can trigger scripts to be written and run from here.

# Import necessary modules
from subprocess import run
from os.path import dirname, abspath, join
from os import makedirs, environ
import sys
from inspect import getfile, currentframe

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir,'src/'))
from IO import getArgs, getFileInfo
import writeProfiles
import writeNamelist
import writeBatch

# Get command line arguments
args = getArgs()

# Organize the directories that we will work in
IOlists = {'profilesIn':args.profilesIn, 'eqIn':args.eqIn, 'saveLoc':args.saveLoc}
maxLen = max([len(data) for key,data in IOlists.items()]) # Due to the checks performed on args, each list will have length maxLen or 1.
longLists = []
for key,data in IOlists.items():
    if len(data) != maxLen:
        IOlists[key] = data * maxLen
    else:
        longLists.append(key)

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
    _, _, _, outDir, _ = getFileInfo('arbitrary/path', actualSaveLoc, 'arbitrary')
    makedirs(outDir, exist_ok=True) # Note that this script has file overwrite powers!

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

    print('***StelloptPlusSfincs: Setup and/or run task(s) {} of {} completed in {}.***'.format(i+1, maxLen, outDir))
