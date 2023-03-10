# This script creates a job.sfincsScan batch script.

def run(profilesInUse, saveLocUse):
    
    '''
    The inputs are set by a wrapper script.
    '''
    
    # Import necessary modules
    from os.path import join
    from os import environ
    from IO import getRunArgs, getFileInfo, writeFile

    # Get command line arguments
    args = getRunArgs()

    # Name output file
    _, _, _, _, outFile = getFileInfo(profilesInUse, saveLocUse, 'job.sfincsScan')

    # Load location of SFINCS directory
    sfincsLoc = join(environ['SFINCS_PATH'],'fortran/version3/sfincs')

    # Load machine name
    machineVar = 'MACHINE'
    machine = environ[machineVar]

    # Create string to be written
    stringToWrite = '#!/bin/bash -l\n'
    
    stringToWrite += '# Standard output and error:\n'
    stringToWrite += '#SBATCH -o ./sfincsJob.out.%j\n'
    stringToWrite += '#SBATCH -e ./sfincsJob.err.%j\n'
    stringToWrite += '# Initial working directory:\n'
    stringToWrite += '#SBATCH -D ./\n'
    stringToWrite += '#\n'
    
    stringToWrite += '# Resource allocation:\n'
    if args.nNodes[0] is not None:
        stringToWrite += '#SBATCH --nodes={}\n'.format(args.nNodes[0])
    if args.nTasksPerNode[0] is not None:
        stringToWrite += '#SBATCH --ntasks-per-node={}\n'.format(args.nTasksPerNode[0])
    if args.nTasks[0] is not None:
        stringToWrite += '#SBATCH --ntasks={}\n'.format(args.nTasks[0])
    if args.mem[0] is not None:
        stringToWrite += '#SBATCH --mem={}\n'.format(args.mem[0])
    stringToWrite += '#\n'
    
    try: # Set up job notification emails if possible 
        email = environ['SFINCS_BATCH_EMAIL']
        
        if args.notifs[0] == 'bad':
            stringToWrite += '#SBATCH --mail-type=fail,invalid_depend,requeue,stage_out\n'
        elif args.notifs[0] == 'all':
            stringToWrite += '#SBATCH --mail-type=all\n'
        elif args.notifs[0] == 'none':
            stringToWrite += '#SBATCH --mail-type=none\n'
        
        stringToWrite += '#SBATCH --mail-user={}\n'.format(email)
        stringToWrite += '#\n'
    
    except KeyError:
        pass
    
    stringToWrite += '# Wall clock limit:\n'
    stringToWrite += '#SBATCH --time={}\n'.format(args.time[0].strip())
    stringToWrite += '\n'
    stringToWrite += '# Load necessary modules (typically must be the same as those used for compiling the code):\n'
    stringToWrite += 'module purge\n'

    if machine == 'raven' or machine == 'cobra':
        stringToWrite += 'module load git\n'
        stringToWrite += 'module load intel/19.1.2\n'
        stringToWrite += 'module load mkl\n'
        stringToWrite += 'module load impi/2019.8\n'
        stringToWrite += 'module load hdf5-mpi/1.10.6\n'
        stringToWrite += 'module load netcdf-mpi/4.7.0\n'
        stringToWrite += 'module load fftw-mpi\n'
        stringToWrite += 'module load anaconda/3/2020.02\n'
        stringToWrite += 'module load petsc-real/3.13.5\n'
        stringToWrite += 'module load mumps-32-noomp/5.1.2\n'
        stringToWrite += '\n'
    
    else:
        from os.path import abspath
        from inspect import getfile, currentframe
        thisFile = abspath(getfile(currentframe()))
        _, thisFileName, _, _, _ = getFileInfo(thisFile, 'arbitrary/path', 'arbitrary')
        errString = 'This machine (as identified by the "{}" environment variable)'.format(machineVar)
        errString += ' is not recognized. Please add the necessary "module load" commands to {}.'.format(thisFileName)
        raise OSError(errString)

    stringToWrite += '# Run the program:\n'
    stringToWrite += 'srun {} -ksp_view\n'.format(sfincsLoc)

    # Write job.sfincsScan file
    writeFile(outFile, stringToWrite)
