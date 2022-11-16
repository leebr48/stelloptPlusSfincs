# This script sets up SFINCS runs and queues them using Slurm

# Import necessary modules
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys
from os import environ
from subprocess import run

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getRunArgs, adjustInputLengths, makeDir, messagePrinter, saveTimeStampFile
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

# Sort out the symmetry argument for *.bc files
if len(args.bcSymmetry) == 1:
    bcSym = args.bcSymmetry * maxLen
else:
    bcSym = args.bcSymmetry

# Loop through the working directories
for i in range(maxLen):
    profilesInUse = IOlists['profilesIn'][i]
    eqInUse = IOlists['eqIn'][i]
    actualSaveLoc = dirname(saveDefaultTarget[i])
    bcSymUse = bcSym[i]
    logString = 'The following automation tasks were carried out:\n'
    appendor = ' file was written\n'

    # Make target directory if it does not exist
    outDir = makeDir(actualSaveLoc) # Note that this script has file overwrite powers!

    # Write requested files
    if not args.noProfiles:
        writeProfiles.run(profilesInUse, outDir)
        logString += '\tprofiles' + appendor

    if not args.noNamelist:
        writeNamelist.run(profilesInUse, outDir, eqInUse, bcSymUse)
        logString += '\tinput.namelist' + appendor

    if not args.noBatch:
        writeBatch.run(profilesInUse, outDir)
        logString += '\tjob.sfincsScan' + appendor
    
    # Call sfincsScan if requested
    if not args.noRun:
        execLoc = join(environ['SFINCS_PATH'],'fortran/version3/utils/sfincsScan')
        cmd = [execLoc]
        userConf = 'with'
        if args.noConfirm:
            cmd.append('arbitraryCommandLineArg')
            userConf = 'without'

        run(cmd, cwd=outDir)
        logString += '\tsfincsScan attempted to run {} user confirmation\n'.format(userConf)
    
    # Save a timestamp file if appropriate
    logString += 'at this time:\n\t'
    saveTimeStampFile(outDir, 'automatedSetupLog', logString)

    # Alert the user that their specified tasks were completed
    messagePrinter('Setup and/or run task(s) {} of {} completed in {}.'.format(i+1, maxLen, outDir))
