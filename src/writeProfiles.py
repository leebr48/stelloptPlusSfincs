# This script creates a SFINCS-readable profiles file.

def run(profilesInUse, saveLocUse):

    '''
    The inputs are set by a wrapper script.
    '''

    # Import necessary modules
    import numpy as np
    from os.path import join
    from matplotlib.pyplot import subplots
    from IO import getRunArgs, getFileInfo, cleanStrings, listifyBEAMS3DFile, makeProfileNames, extractProfileData, sortProfileFunctions, generatePreamble, generateDataText, writeFile, messagePrinter, prettyRadialVar
    from dataProc import findMinMax, scaleInputData, nonlinearInterp

    # Get command line arguments
    args = getRunArgs()

    # Name input and output files
    inFile, _, _, outDir, outFile = getFileInfo(profilesInUse, saveLocUse, 'profiles') # Name mandated by SFINCS
    
    plotName = 'interpFuncFit'
    plotFile = join(outDir, plotName+'.pdf')

    # Clean input variable names and do some clerical checks
    varsToFind = ['NE', 'NI', 'TE', 'TI']
    if args.loadPot:
        varsToFind.append('POT')
    
    prefixesOfInterest = cleanStrings(varsToFind)

    # Extract the data from the BEAMS3D input file
    listifiedInFile = listifyBEAMS3DFile(inFile)

    varsOfInterest = makeProfileNames(prefixesOfInterest)
    dataOfInterest = extractProfileData(listifiedInFile, varsOfInterest)

    radialBounds = findMinMax(dataOfInterest)
    
    # Scale the data according to the reference variable values
    scaledData = scaleInputData(dataOfInterest)

    # Interpolate the data in case the radial lists do not all contain the same points
    ders = {}
    for key,val in scaledData.items():
        ders[key] = 0

    if args.loadPot:
        ders['pot'] = 1 # Only take a derivative when we'll need it for further calculations

    interpolatedData = nonlinearInterp(scaledData, ders, k=3)
    sortedInterpolatedData = sortProfileFunctions(interpolatedData) # Note that "pot" is not included in this output even if it is included in the input

    # Gather the components of profiles file
    radial_coordinate_ID = 1 # Corresponds to normalized toroidal flux, which is "s" in STELLOPT and "psiN" in SFINCS

    radii = list(np.linspace(start=radialBounds['min'], stop=radialBounds['max'], num=args.numInterpSurf[0], endpoint=True))

    # Note that NErs, generalEr_min, and generalEr_max are only used by SFINCS if scanType = 5.
    NErs = lambda x: args.numErSubscan[0]
    
    if args.loadPot:
        generalEr_min = lambda x: interpolatedData['pot'][0](x) + args.minSeedEr[0]
        generalEr_max = lambda x: interpolatedData['pot'][0](x) + args.maxSeedEr[0]
    else:
        generalEr_min = lambda x: args.minSeedEr[0]
        generalEr_max = lambda x: args.maxSeedEr[0]

    funcs = [NErs, generalEr_min, generalEr_max]
    funcs.extend(sortedInterpolatedData)

    # Plot the fitted interpolation functions to ensure they represent the data well
    fig,ax = subplots()

    if args.loadPot:
        scaledData.pop('pot') # We don't have data for the *derivative* of Phi, so this graph looks out of place

    leg = []
    for key, data in scaledData.items():

        for specInd, (IVvec, DVvec) in enumerate(zip(data['iv'], data['dv'])):
            color = next(ax._get_lines.prop_cycler)['color']
            ax.scatter(IVvec, DVvec, c=color)
            ax.plot(radii, interpolatedData[key][specInd](radii), c=color)
            keyUse = key + str(specInd+1)
            leg.append(keyUse)

    ax.legend(leg, loc='best')
    ax.set_xlabel(prettyRadialVar('psiN'))
    ax.set_ylabel('Normalized Value')

    fig.savefig(plotFile, bbox_inches='tight', dpi=400)
    messagePrinter('{} plot created.'.format(plotName))

    # Get the string to write in profiles file
    stringToWrite = generatePreamble(radial_coordinate_ID)
    stringToWrite += generateDataText(radii, *funcs)

    # Write profiles file
    writeFile(outFile, stringToWrite)
