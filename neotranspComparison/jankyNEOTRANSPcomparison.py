# NOTE: the electron and ion ordering in all the loaded files must be the same for this script to work!!!

# User inputs
baseDirs = ['/u/lebra/src/stelloptPlusSfincs/neotranspComparison/sfincs']
plotTitles = False # Can be True, False, or a string

# Imports
from os.path import basename, join
import numpy as np
import matplotlib.pyplot as plt

# Useful functions
def getData(baseDir, sfincsFileExt, neotranspFileExt, compareFileExt):
    
    runBaseName = basename(baseDir)
    sfincsDir = join(baseDir, 'processed')
    neotranspDir = join(baseDir, 'neotransp')
    saveDir = neotranspDir

    sfincsFile = join(sfincsDir, runBaseName + '-{}-vs-rN.dat'.format(sfincsFileExt))
    neotranspFile = join(neotranspDir, 'input.' + runBaseName + '.{}_vs_rho'.format(neotranspFileExt))

    figFile = join(saveDir, runBaseName + '_{}_COMPARE.png'.format(compareFileExt))
    
    sfincsData = np.loadtxt(sfincsFile)
    neotranspData = np.loadtxt(neotranspFile)

    return runBaseName, sfincsDir, figFile, sfincsData, neotranspData

def makeMainPlot(sfincsData, neotranspData, neotranspDataMultip, ylabel, title=None):
    
    plt.figure()

    plt.plot(sfincsData[:,0], sfincsData[:,1:])
    plt.plot(neotranspData[:,0], neotranspData[:,1:] * neotranspDataMultip)

    if title is not None:
        plt.title(title)
    plt.xlabel(r'$r_{N} = \rho = r_{\mathrm{eff}}/a_{\mathrm{eff}}$')
    plt.ylabel(ylabel)

def makeLegendAndSave(sfincsData, sfincsDir, runBaseName, sfincsFileExt, figFile):

    defaultleg = ['sfincs', 'neotransp']
    if sfincsData.shape[1] > 2:
        ZsFile = join(sfincsDir, runBaseName + '-{}-vs-rN.Zs'.format(sfincsFileExt))
        Zs = np.loadtxt(ZsFile)
        leg = []
        for software in defaultleg:
            for Z in Zs:
                leg.append(software + ', ' + r'$Z={}$'.format(int(Z)))
    else:
        leg = defaultleg
    
    plt.legend(leg, loc='best')

    plt.savefig(figFile, bbox_inches='tight', dpi=400)
    plt.close('all')

def makePlot(baseDir, sfincsFileExt, neotranspFileExt, neotranspDataMultip, vertAxisLabel, compareFileExt, title=True):
    runBaseName, sfincsDir, figFile, sfincsData, neotranspData = getData(baseDir, sfincsFileExt, neotranspFileExt, compareFileExt)
    if title:
        if isinstance(title,str):
            titleToUse = title
        else:
            titleToUse = runBaseName
    else:
        titleToUse = None
    makeMainPlot(sfincsData, neotranspData, neotranspDataMultip, vertAxisLabel, titleToUse)
    makeLegendAndSave(sfincsData, sfincsDir, runBaseName, sfincsFileExt, figFile)

# Actually make the plots
for baseDir in baseDirs:

    # Radial electric field
    makePlot(baseDir, 'Er', 'ErkVm', 1000, r'Radial electric field $\mathrm{\left(\frac{V}{m}\right)}$', 'Er', title=plotTitles)

    # Neoclassical particle flux
    makePlot(baseDir, 'extensiveNeoclassicalParticleFlux', 'ParticlefluxEI', 1, r'Neoclassical particle flux $\mathrm{\left(\frac{1}{s}\right)}$', 'NeoclassicalParticleFlux', title=plotTitles)

    # Neoclassical heat flux
    makePlot(baseDir, 'extensiveNeoclassicalHeatFlux', 'HeatfluxEI', 10**6, r'Neoclassical heat flux $\mathrm{\left(\frac{J}{s}\right)}$', 'NeoclassicalHeatFlux', title=plotTitles)

    # Bootstrap current
    makePlot(baseDir, 'FSABjHatOverB0', 'Bootstrapcurrdens', 1, r'Bootstrap current $\mathrm{\left(\frac{A}{m^{2}}\right)}$', 'BootstrapCurrent', title=plotTitles)