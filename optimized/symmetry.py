from configuration import Configuration
from typing import List
from rotation import rotate_configuration
import numpy as np
from utils import write_to_json
from distance import calculate_distance
from displacement import AXES, Displacement


def create_configuration_dict(configurations:List[Configuration]):
  result = {}
  for configuration in configurations:
    result[configuration.key] = configuration
  return result


def find_symmetries_for_configuration(configuration, configurations, initial_symmetries, operators, tol):
  configurations_dict = create_configuration_dict(configurations)
  symmetries = {}

  for initial_sym_atom_index, operators_indices in initial_symmetries.items():
    for operator_index in operators_indices:
      operator = operators[operator_index]
      rotated_configuration = rotate_configuration(configuration, operator) # Obviously don't need to repeat this keep somewhere
      rotated_scaled_position = rotated_configuration.atoms.get_scaled_positions()[configuration.atom_index]
      
      for axis in AXES:
        for sgn in ["-", "+"]:
           alt_configuration_key = Configuration.create_key(initial_sym_atom_index, axis, sgn)
           alt_configuration = configurations_dict[alt_configuration_key]
           alt_scaled_position = alt_configuration.atoms.get_scaled_positions()[alt_configuration.atom_index]
           distance = calculate_distance(rotated_scaled_position, alt_scaled_position)
           if distance <= tol:
             if configuration.key in symmetries:
               if alt_configuration.key in symmetries[configuration.key]:
                 symmetries[configuration.key][alt_configuration.key].append(operator_index)
               else:
                 symmetries[configuration.key][alt_configuration.key] = [operator_index]
             else:
               symmetries[configuration.key] = {}
               symmetries[configuration.key][alt_configuration.key] = [operator_index]       
  return symmetries

def find_symetries(configurations:List[Configuration], operators, mapping, tol = 1e-4):
  symmetries = {}
  for i, configuration in enumerate(configurations):
    if configuration.atom_index in mapping:
      initial_symmetries = mapping[configuration.atom_index]
      preserved_symmetries = find_symmetries_for_configuration(configuration, configurations, initial_symmetries, operators, tol)
      symmetries[configuration.key] = preserved_symmetries[configuration.key]
  write_to_json(symmetries, './data/mapping/displaced_atom_mappings.json')
  return symmetries


def symmetrize_equilibrium_potential(eq_potential, eq_potential_mapping):
  v = eq_potential['data'][0].ravel()
  v_sym = np.zeros(v.shape)
  for operator in eq_potential_mapping:
    v_sym += v[eq_potential_mapping[operator]] 
  v_sym = v_sym/len(eq_potential_mapping.keys())
  return {'displacement_key': 'eq' , 'data' : v_sym.reshape(eq_potential['data'].shape)}


def symmetrize_displaced_potential(displaced_potential, eq_potential_mapping, displaced_atom_mapping):
  displacement_key = displaced_potential['displacement_key']
  shape = displaced_potential['data'].shape
  v = displaced_potential['data'][0].ravel()
  v_sym = np.zeros(v.shape)
  if displacement_key in displaced_atom_mapping:
    symmetries = displaced_atom_mapping[displacement_key]
    for atom in symmetries:
      operators = symmetries[atom]
      for operator in operators:
        v_sym += v[eq_potential_mapping[operator]]
      v_sym = v_sym/len(operators)
    v_sym = v_sym/len(symmetries.keys())

  return {"displacement_key": displacement_key, "data": v_sym.reshape(shape)}