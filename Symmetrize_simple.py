### after the first step symmetrize the slab using the following code 
## this might need some changing depending on the system 
#!/usr/bin/env python3
import numpy as np
import ase.io
import ase.constraints
from copy import deepcopy
from ase import Atoms,Atom
atoms = ase.io.read('CONTCAR')
atoms.center(axis=2)
ase.io.write('take.traj', atoms)
base_z = np.array([atom.z for atom in atoms]).min()
base_layer_index = [atom.index for atom in atoms if abs(base_z - atom.z) < 0.1]
base_layer_center = np.array([atom.position for atom in atoms[base_layer_index]]).mean(axis=0)
print(base_layer_center)
inverted_atoms = Atoms([deepcopy(atom) for atom in atoms if atom.index not in base_layer_index])
for atom in inverted_atoms:
    atom.position = 2 * base_layer_center - atom.position
atoms += inverted_atoms
del atoms.constraints
constraints = ase.constraints.FixAtoms(indices=base_layer_index)
atoms.set_constraint(constraints)
cell = atoms.get_cell()
cell[2][2] = 60
atoms.set_cell(cell)
atoms.center(axis=2)
ase.io.write('sym.traj', atoms)
