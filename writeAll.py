# This script writes both the profiles and input.namelist files.
# It is simply a wrapper for the scripts that perform these tasks individually.

from IO import getArgs

# Get command line arguments
args = getArgs()

if not args.noProfiles:
    with open('writeProfiles.py', mode='r', encoding='utf-8') as code:
        exec(code.read())

if not args.noNamelist:
    with open('writeNamelist.py', mode='r', encoding='utf-8') as code:
        exec(code.read())
