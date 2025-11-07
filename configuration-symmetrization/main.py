import numpy as np
from gpaw import GPAW
from displacement import get_all_displacements
from configuration import create_configurations
from symmetry import find_symetries
from mapping import create_equilibrium_atom_mapping_table

delta = 1e-3
displacements = get_all_displacements(delta)

cell_type = "supercell"
calc = GPAW(f'./data/si_{cell_type}.gpw', parallel={'domain': 1, 'band': 1})
operators = np.load(f'./data/si_{cell_type}_symmetries.npy')
atoms = calc.atoms
equilibrium_symmetries = create_equilibrium_atom_mapping_table(atoms, operators)


configurations = create_configurations(atoms, displacements)
find_symetries(configurations, operators, equilibrium_symmetries)
# print(len(configurations))