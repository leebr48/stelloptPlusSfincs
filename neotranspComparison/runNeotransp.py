# This script runs Neotransp, but it does NOT make the BEAMS3D inputs
# For now, the relevant parts will be as close to the Matlab script as possible.
# FIXME maybe make this work with impurities?

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

lib_path=os.getenv('NEOTRANSP_PYTHON_LIB')
sys.path.insert(0,lib_path)
from W7Xdatabase import neoclass_database
from neolib import profile_data, transp_data

# Can change
ldeuterium = False # Inject D into plasma or not
ne_vec = [5E19] # m^-3
te_vec = [1000] # eV
conf_name='EIM'; equilID = 'w7x-sc1'
a = 0.51092 # Effective minor radius of W7-X
B00axis = 2.41 # Hack according to Hakan Smith... you could also do Baxis_phi0=2.52

# Were hard-coded into Matlab, probably shouldn't change
ti_max = 1.6E3 # eV
B = 2.5

#####################

h0 = 0.05
rho = np.linspace(start=h0, stop=1-h0, num=int(1/h0-1))
s = rho ** 2
r = a * np.sqrt(s) # As defined in the Neotransp docs

# Profiles
f_ne = lambda s: ne_max / B * (B - s**3)
f_te = lambda s: te_max * (1 - s**0.5)
f_ti = lambda s: np.where(f_te(s) <= ti_max, f_te(s), ti_max)

# Profile derivatives
dfds_ne = lambda s: -ne_max / B * (3 * s**2)
dfds_te = lambda s: -0.5 * te_max * s**(-0.5)
dfds_ti = lambda s: np.where(f_te(s) <= ti_max, dfds_te(s), 0) # This was a numerical derivative in Matlab

# Converting to derivative wrt r for Python version of Neotransp
dsdr = lambda s: 2 * np.sqrt(s) / a
dfdr_ne = lambda s: dfds_ne(s) * dsdr(s)
dfdr_te = lambda s: dfds_te(s) * dsdr(s)
dfdr_ti = lambda s: dfds_ti(s) * dsdr(s)

for (ne_max, te_max) in zip(ne_vec, te_vec):
    
    # File name
    ext = 'W7X_' + conf_name + '_n' + '%2.2i'%(ne_max/1E19) + '_e' + '%2.2i'%(te_max/100) + '_i' + '%2.2i'%(f_ti(0)/100)
    if ldeuterium:
        ext += '_8SD2'
    else:
        ext += '_8SH2'
    
    # Load profiles into Neotransp
    Prof = profile_data(['e', 'H'], r=r) # Class initialization
    Prof.set_n_20m3('e', f_ne(s) / 10**20)
    Prof.set_T_keV('e', f_te(s) / 10**3)
    Prof.set_dn_20m3dr('e', dfdr_ne(s) / 10**20)
    Prof.set_dT_keVdr('e', dfdr_te(s) / 10**3)
    Prof.set_n_20m3('H', f_ne(s) / 10**20) # Note that ne = ni
    Prof.set_T_keV('H', f_ti(s) / 10**3)
    Prof.set_dn_20m3dr('H', dfdr_ne(s) / 10**20) # Note that ne = ni
    Prof.set_dT_keVdr('H', dfdr_ti(s) / 10**3)

    # Get DKES and Boozer data
    db = neoclass_database(lib_path+'/w7xOp1.2.nc') # Can change the .nc file to "w7x.nc" if desired... superset of current data, but very slow
    dk = db.extract(equilID)

    # Check derivatives just to be safe
    absthreshold = 1E-3
    gradCheck = Prof.validate_derivatives(dn_20m3dr_absthreshold=absthreshold, dT_keVdr_absthreshold=absthreshold)
    profilesCheckFileName = 'profilesCheck.png'
    fig,ax = Prof.plot(withgradients='lin', withcollisionality=True, dk=dk, additional_prof='grad approx', savefile='./{}'.format(profilesCheckFileName))

    if gradCheck is False:
        stringToPrint = 'At least one of the input gradients does not agree with the numerical gradient derived from cubic spline interpolation (absolute difference > {}).'.format(absthreshold)
        stringToPrint += 'This could be an error or a numerical artifact. Check {} to see.'.format(profilesCheckFileName)
        print(stringToPrint)

    # Run Neotransp
    Transp = transp_data(Prof, dk=dk, B00=B00axis, roots='i/e')
    
    # Save outputs
    transpOutputFileName = 'transpResults'
    Transp.plot(totalflux=True, xlabel='r/a', showsum_bootstrap=True, savefile=transpOutputFileName+'.pdf')
    Transp.makeSFINCSscan21runspec('runspec.dat')
    
    #Transp.save(transpOutputFileName+'.txt', quantities=['rho', 'Er', 'Jbs', 'Flux', 'EnergyFlux'], fluxmode='total', energyfluxmode='total')
    toSave = np.c_[rho, Transp.get_ErkVm()]
    np.savetxt('input.ErkVm_vs_rho', toSave)
    toSave = np.c_[rho, Transp.get_Jbs_kAm2()]
    np.savetxt('input.Bootstrapcurrdens_vs_rho', toSave)
    toSave = np.c_[rho, Transp.get_NCflux('e', mode='total'), Transp.get_NCflux('H', mode='total')]
    np.savetxt('input.ParticlefluxEI_vs_rho', toSave)
    toSave = np.c_[rho, Transp.get_NCenergyflux('e', mode='total'), Transp.get_NCenergyflux('H', mode='total')]
    np.savetxt('input.HeatfluxEI_vs_rho', toSave)
