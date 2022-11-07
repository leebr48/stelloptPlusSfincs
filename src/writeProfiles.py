# This script creates a SFINCS-readable profiles file.

def run(profilesInUse, saveLocUse):

    '''
    The inputs are set by a wrapper script.
    '''

    # Import necessary modules
    # FIXME double-check at the end to ensure all of these imports are still needed!
    import numpy as np
    from os.path import join
    from matplotlib.pyplot import subplots
    from IO import getRunArgs, getFileInfo, cleanStrings, listifyBEAMS3DFile, extractProfileData, makeProfileNames, generatePreamble, generateDataText, writeFile, messagePrinter
    from dataProc import findMinMax, scaleInputData, nonlinearInterp

    # Get command line arguments
    args = getRunArgs()

    # Name input and output files
    inFile, _, _, outDir, outFile = getFileInfo(profilesInUse, saveLocUse, 'profiles') # Name mandated by SFINCS
    
    plotName = 'interpFuncFit'
    plotFile = join(outDir, plotName+'.pdf')

    # Clean input variable names and do some clerical checks
    prefixesOfInterest = cleanStrings(['NE', 'NI', 'TE', 'TI'])

    # Extract the data from the BEAMS3D input file
    listifiedInFile = listifyBEAMS3DFile(inFile)

    varsOfInterest = makeProfileNames(prefixesOfInterest)
    dataOfInterest = extractProfileData(listifiedInFile, varsOfInterest)
    print(dataOfInterest); print('**********************')

    radialBounds = findMinMax(dataOfInterest)
    print(radialBounds); print('**********************')
    
    # Scale the data according to the reference variable values
    scaledData = scaleInputData(dataOfInterest)
    print(scaledData); print('**********************')
    # Interpolate the data in case the radial lists do not all contain the same points
    interpolatedData = nonlinearInterp(scaledData) # FIXME at some point, you will need to set the temperature profiles to be the same for each ion if they're not specified... maybe even before this, in a different function?
    print(interpolatedData); quit()

    # Gather the components of profiles file
    radial_coordinate_ID = 1 # Corresponds to normalized toroidal flux, which is S in STELLOPT and psiN in SFINCS

    radii = list(np.linspace(start=radialBounds['min'], stop=radialBounds['max'], num=args.numInterpSurf[0], endpoint=True))

    # Note that NErs, generalEr_min, and generalEr_max are only used by SFINCS if scanType = 5.
    NErs = lambda x: args.numManErScan[0]
    generalEr_min = lambda x: args.minEr[0]
    generalEr_max = lambda x: args.maxEr[0]

    funcs = [NErs, generalEr_min, generalEr_max]
    # FIXME everything works up to here, I think
    funcs.extend([interpolatedData[prefix] for prefix in prefixesOfInterest]) #FIXME you need to unpack interpolatedData properly

    # Plot the fitted interpolation functions to ensure they represent the data well
    fig,ax = subplots()

    leg = []
    for key, data in scaledData.items():
        ax.scatter(data['iv'], data['dv'])
        ax.plot(radii, interpolatedData[key](radii)) #FIXME interpolatedData probably won't work like this anymore
        leg.append(key)

    ax.legend(leg, loc='best')
    ax.set_xlabel(r'SFINCS $\psi_{N}$ $\left(= \mathrm{STELLOPT}{\ }S\right)$')
    ax.set_ylabel('Normalized Value')

    fig.savefig(plotFile, bbox_inches='tight', dpi=400) # FIXME ensure that the plot is still correct (that is, your interpolation functions are still working)
    messagePrinter('{} plot created.'.format(plotName))

    # Get the string to write in profiles file
    stringToWrite = generatePreamble(radial_coordinate_ID)
    stringToWrite += generateDataText(radii, *funcs)

    # Write profiles file
    writeFile(outFile, stringToWrite)
