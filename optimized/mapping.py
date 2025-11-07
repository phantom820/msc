import numpy as np
from rotation import rotate_atom, rotate_scaled_positions
from distance import calculate_distance
from utils import scaled_positions_to_grid_positions, grid_to_scaled_positions, write_to_json

def find_atom_closest_match(index, rotated_atoms, atoms, epsilon = 1e-3):
  rotated_scaled_position = rotated_atoms.get_scaled_positions()[index]
  ref_positions = atoms.get_scaled_positions()
  distance = calculate_distance(rotated_scaled_position, ref_positions, epsilon)
  idx = int(np.argmin(distance))
  return idx, distance[idx]

def add_equilibrium_atom_mapping(mapping, atom_index, alt_atom_index,  operator_index):
  if atom_index in mapping:
    if alt_atom_index in mapping[atom_index]:
      mapping[atom_index][alt_atom_index].append(operator_index)
    else:
      mapping[atom_index][alt_atom_index] = [operator_index]
  else:
    mapping[atom_index] = {}
    mapping[atom_index][alt_atom_index] = [operator_index]

def verify_equilibrium_potential_mapping(mapping_idx):
  unique, counts = np.unique(list(mapping_idx), return_counts = True)
  assert len(mapping_idx) == unique.shape[0], "Not all indices appear in mapping"
  assert np.all(counts == 1) == True, "Duplicate index exists in mapping"

def verify_equilibrium_atom_mapping(mapping, number_of_operators):
  for _ in mapping.keys():
    temp = 0
    for k,v in mapping[_].items():
      temp = temp + len(v)
    
    assert temp == number_of_operators, "Equillibrium system missing mapping for some operators"

def create_equilibrium_atom_mapping_table(atoms, rotations, tol = 1e-4):
  mapping = {}
  for i in range(len(atoms)):
    for j, R in enumerate(rotations):
      rotated_atoms = rotate_atom(i, atoms, R)
      idx, distance = find_atom_closest_match(i, rotated_atoms, atoms)
      if distance <= tol:
        add_equilibrium_atom_mapping(mapping, i, idx, j)
  verify_equilibrium_atom_mapping(mapping, len(rotations))
  write_to_json(mapping, "./data/mapping/equilibrium_atom_mappings.json")
  return mapping

def create_equilibrium_potential_mapping_table(rotations, shape):
  mapping = {}
  scaled_positions = grid_to_scaled_positions(shape)
  for i in range(len(rotations)):
    R = rotations[i]
    rotated_scaled_positions = rotate_scaled_positions(scaled_positions, R)
    multi_idx = scaled_positions_to_grid_positions(rotated_scaled_positions, shape)
    flat_idx = np.ravel_multi_index(multi_idx.T, shape).astype(int)
    verify_equilibrium_potential_mapping(flat_idx)
    mapping[i] = flat_idx.tolist()
  write_to_json(mapping, "./data/mapping/equilibrium_potential_mappings.json")
  return mapping