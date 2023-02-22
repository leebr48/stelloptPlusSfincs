# FIXME this simple script computes the collisionality of particles at given conditions

from scipy.constants import m_e # Electron mass in kg
from scipy.constants import m_p # Proton mass in kg
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from dataProc import nu_ab, approx_nu_ii

# Sort out inputs
ne = 1.8 # FIXME make arg
te = 11 # FIXME make arg
nis = [0.8, 0.09] # FIXME make arg
tis = [11, 11] #FIXME make arg
zis = [2, 2] # FIXME make arg
mis = [2, 3] # FIXME make arg
K = 1 # FIXME make arg

# Set up other quantities
ze = -1
me = m_e / m_p

# Set up electron dictionary
eDict = {'n':ne, 't':te}

# Perform nu_ee calculations first # FIXME actually probably only calculate ion stuff?... Your functions might not work for nu_ei (vs nu_ie)
aDict = {'n':ne, 't':te, 'z':ze, 'm':me}
bDict = {'n':ne, 't':te, 'z':ze, 'm':me}

nu_ee = nu_ab(aDict, bDict, K)
print(nu_ee)

# FIXME your approximation function for nu may be broken, or it may only work for electron-ion versus the other way?
# FIXME you can probably fix your coulomb logarithm stuff based on neotransp, and just make everything 'species a and b'
# FIXME if you are inputing ion species twice, probably check that their inputs are truly identical (or just force it) so the Coulomb log calc is correct

#for nb, tb, zb, mb in zip(nis, tis, zis, mis):
