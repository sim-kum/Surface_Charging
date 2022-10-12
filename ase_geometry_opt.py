##### first step is to do local geometry optimization without solvation or surface charging ### 

from ase.constraints import Hookean
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

os.environ["VASP_PP_PATH"]="/home/ksimran/PP/pp_new/POT/"

try:
    with open('nelect.pickle', 'rb') as picklefile:
        nelect_to_update = pickle.load(picklefile)
except:
    nelect_to_update = {}

try:
    slab = ase.io.read('opt.traj')
except:
    try:
        slab = ase.io.read('input.traj')
    except:
        slab = ase.io.read('CONTCAR')
calc_1 = Vasp(xc='pbe',
         ldau_luj={'Co':{'L' : 2, 'U': 3.52,'J' :0}},  ##specify any U parameters
         encut=400, ## lower encut values
         kpts=[1,1,1], ## local optimization with only gama point
         ediff = 1.00e-06,
         ediffg = -0.1, ## very low convergance
         prec='A',
         algo='fast',
         ispin=1, #non spin polarized calculation
         ismear=2, #using the second order MP smearing
         sigma=0.1,
         nelm = 60,
         nsw = 300,
         ibrion = 2,
         isif = 2,
         lwave = False,
         lcharg = False,
         npar = 4,
         lreal='A', #projection done in real space
         kgamma=True,
         istart=0, #if you want to restart put other values
         icharg=2)
slab.set_calculator(calc_3)
slab.get_potential_energy()
##Repeat calculation with stringent convergance criteria
##Usually takes a few steps
slab = ase.io.read('CONTCAR')
calc_2 = Vasp(xc='pbe',
         ldau_luj={'Co':{'L' : 2, 'U': 3.52,'J' :0}},
         encut=400,
         kpts=[6,3,1],
         ediff = 1.00e-06,
         ediffg = -0.05,
         prec='A',
         algo='fast',
         ispin=2, #non spin polarized calculation
         ismear=2, #using the second order MP smearing
         sigma=0.1,
         nelm = 60,
         nsw = 300,
         ibrion = 2,
         isif = 2,
         lwave = False,
         lcharg = False,
         npar = 4,
         lreal='A', #projection done in real space
         kgamma=True,
         istart=0, #if you want to restart put other values
         icharg=2)
slab.set_calculator(calc_2)
slab.get_potential_energy()
slab.write('optimized.traj')
