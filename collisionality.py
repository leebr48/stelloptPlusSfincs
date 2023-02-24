# This simple script computes the collisionality of the input particles. It is assumed that all the particles live together in the same plasma.

# Load necessary modules
from scipy.constants import m_e # Electron mass in kg
from scipy.constants import m_p # Proton mass in kg
from os.path import dirname, abspath, join
from inspect import getfile, currentframe
import sys

thisDir = dirname(abspath(getfile(currentframe())))
sys.path.append(join(thisDir, 'src/'))
from IO import getCollisionalityArgs, messagePrinter
from dataProc import nu_ab, thermalVelocity

# Sort out inputs
args = getCollisionalityArgs()
ns = args.ns
ts = args.ts
zs = args.zs
ms = [m_e / m_p if x < 0 else x for x in args.ms]
Ks = args.Ks

# Perform calculations
nus = []
nu_vs = []
for na, ta, za, ma, Ka in zip(ns, ts, zs, ms, Ks):
    aDict = {'n':na, 't':ta, 'z':za, 'm':ma}
    nua = 0
    for nb, tb, zb, mb in zip(ns, ts, zs, ms):
        bDict = {'n':nb, 't':tb, 'z':zb, 'm':mb}
        nuab = nu_ab(aDict, bDict, Ka)
        nua += nuab
    nua_va = nua / thermalVelocity(ta * Ka, ma, units='keV_mp')
    nus.append(nua)
    nu_vs.append(nua_va)

# Report results
nus_msg = 'The collisionalities (nu, in Hz) of the input particles, in order, are:\n'
nus_msg += str(nus)
messagePrinter(nus_msg)

nu_vs_msg = 'The collisionalities (nu/v, in 1/m) of the input particles, in order, are:\n'
nu_vs_msg += str(nu_vs)
messagePrinter(nu_vs_msg)
