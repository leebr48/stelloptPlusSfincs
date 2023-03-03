# stelloptPlusSfincs

Python scripts for making certain aspects of [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT) and [SFINCS](https://github.com/landreman/sfincs) work nicely together. The current focus is on running many SFINCS cases in a consistent and easy way.

## Installation

In principle, the scripts in this repository can perform their core functions given only the appropriate input files. In practice, it is nearly essential that you install SFINCS (instructions in its repository), and it is recommended that you install STELLOPT (instructions [here](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation)) as well. This will make everything work smoothly. SFINCS is memory-intensive, so the use of a high-performance cluster will likely be required unless you only run very small test problems.

The Python 3 [standard library](https://docs.python.org/3/library/index.html) should be installed, as several of its packages are utilized. The [NumPy](https://numpy.org/), [SciPy](https://scipy.org/), [h5py](https://www.h5py.org/), and [Matplotlib](https://matplotlib.org/) packages are also required.

Set the environment variable (preferably permanently, such as in your `.bashrc` file) `SFINCS_PATH=/path/to/sfincs/repository`. If you want to receive job updates from Slurm, also set `SFINCS_BATCH_EMAIL=your_email@website.com` in the same way.

## Use

Currently, these scripts can take BEAMS3D input files and a number of command line arguments and use them to create the files needed for SFINCS to run. SFINCS can also be run automatically, and its outputs can be processed easily. Other scripts are included for convenience. You can see how to use this repository by running any of the scripts in the main directory with the `--help` flag. The scripts themselves also contain notes on their use.

Note that the profiles input into these scripts are not always checked for physical sensibility. They must satisfy quasineutrality, for instance, or the results may not be reliable. In general, the density, temperature, and radial gradients of these quantities must be specified for all species (electrons and ions) on every flux surface for which SFINCS will perform calculations. It is easiest to specify profiles thoughout the plasma volume and let the software calculate the necessary values from them. If desired, you may specify a single electron temperature profile and a single ion temperature profile; the ion temperature profile will be used for all ion species in this case. The masses and charges of all ions must be provided in the standard BEAMS3D format.

## Version Notes

The scripts in this repository have been tested on the `master` branch of SFINCS. At the time of writing, this branch was last updated on 12 May 2022. The version of STELLOPT used is likely not critical since its basic calls and output structure (especially for VMEC) have not changed much over time.
