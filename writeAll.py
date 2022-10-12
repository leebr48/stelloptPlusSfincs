# This script writes files necessary for SFINCS and sfincsScan runs.
# It is simply a wrapper for the scripts that write each file individually.
# Note that exec() commands are used in this script so that the called scripts 
# can still be run individually, if desired.

from IO import getArgs

# Get command line arguments
args = getArgs()

if not args.noProfiles:
    with open('writeProfiles.py', mode='r', encoding='utf-8') as code:
        exec(code.read())

if not args.noNamelist:
    with open('writeNamelist.py', mode='r', encoding='utf-8') as code:
        exec(code.read())

if not args.noBatch:
    with open('writeBatch.py', mode='r', encoding='utf-8') as code:
        exec(code.read())
