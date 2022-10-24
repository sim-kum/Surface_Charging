from ase.build import surface
from ase import Atoms
from ase.build import molecule, add_adsorbate
from ase.optimize import BFGS
from ase.calculators.vasp import Vasp
import os
import subprocess
from ase.calculators.vasp import Vasp2
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.lattice.hexagonal import *
import ase.io
import pickle
try:
    with open('nelect.pickle', 'rb') as picklefile:
        nelect_to_update = pickle.load(picklefile)
except:
    nelect_to_update = {}

try:
    slab = ase.io.read('opt.traj')
except:
    try:
        slab = ase.io.read('POSCAR')
    except:
        slab = ase.io.read('init.traj')

# Make a test slab
calc = Vasp2(xc='pbe',
         encut=400, ## depending on system
         kpts=[5,5,1], ## depending on system 
         ediff = 1.00e-06,
         ediffg = -0.01,
         prec='A',
         algo='fast',
         ispin = 1,
         ismear=2, #using the second order MP smearing
         sigma=0.1, ## depending on system
         nelm = 60,
         nsw = 500,
         ibrion = 2,
         isif = 2,
         npar = 4,
         ivdw = 4,
         lsol = True,tau=0.0,lambda_d_k=9.6,eb_k=78.4,
         lreal='A', #projection done in real space
         kgamma=True,
         istart=1, #if you want to restart put other values
         icharg=2,
         **nelect_to_update)
slab.set_calculator(calc)
slab.get_potential_energy()
