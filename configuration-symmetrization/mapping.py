import numpy as np
from gpaw import GPAW
from gpaw.symmetry import Symmetry
from utils import write_to_json
from rotation import rotate_atom
from distance import calculate_distance

def find_atom_closest_match(index, rotated_atoms, atoms, epsilon = 1e-3):
  rotated_scaled_position = rotated_atoms.get_scaled_positions()[index]
  ref_positions = atoms.get_scaled_positions()
  distance = calculate_distance(rotated_scaled_position, ref_positions, epsilon)
  idx = int(np.argmin(distance))
  return idx, distance[idx]


def add_equilibrium_mapping(mapping, atom_index, operator_index, alt_atom_index):
  if atom_index in mapping:
    mapping[atom_index][operator_index] = alt_atom_index
  else:
    mapping[atom_index] = {}
    mapping[atom_index][operator_index] = alt_atom_index

def create_equilibrium_atom_mapping_table(atoms, operators, tol = 1e-3):
  mapping = {}
  for i, atom in enumerate(atoms):
    for j, R in enumerate(operators):
      rotated_atoms = rotate_atom(i, atoms, R)
      idx, distance = find_atom_closest_match(i, rotated_atoms, atoms)
      if distance <= tol:
        print(distance)
        add_equilibrium_mapping(mapping, i, j, idx)
  write_to_json(mapping,  './data/mapping/equilibrium_atom_mappings.json')
  return mapping

if __name__ == "__main__":
  cell_type = "supercell"
  calc = GPAW(f'./data/si_{cell_type}.gpw', parallel={'domain': 1, 'band': 1})
  symmetries = np.load(f'./data/si_{cell_type}_symmetries.npy')
  atoms = calc.atoms
  rotated_atoms = rotate_atom(2, atoms, symmetries[4])
  # print(atoms.get_positions()[1])
  # print(rotated_atoms.get_positions()[1])
  r = create_equilibrium_atom_mapping_table(atoms, symmetries)
  # create_displaced_atom_mapping_table(atoms, symmetries)
  # print(r)