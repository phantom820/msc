from configuration import Configuration
from typing import List
from rotation import rotate_configuration
import numpy as np
from utils import write_to_json
from displacement import DIRECTIONS, SIGNS


def create_configuration_dict(configurations):
  result = {}
  for configuration in configurations:
    sign = "+" if configuration.delta > 0 else "-"
    if configuration.atom_index in result:
      if configuration.direction in result[configuration.atom_index]:
        result[configuration.atom_index][configuration.direction][sign] = configuration
      else:
         result[configuration.atom_index][configuration.direction] = {}
         result[configuration.atom_index][configuration.direction][sign] = configuration
    else:
       result[configuration.atom_index] = {}
       result[configuration.atom_index][configuration.direction] = {}
       result[configuration.atom_index][configuration.direction][sign] = configuration

  return result


def find_symmetries_for_configuration(configuration, configurations, initial_symmetries, operators, tol):
  initial_symmetries_atom_indices = set(initial_symmetries.values())
  initial_symmetries_operator_indices = [int(i) for i in list(initial_symmetries.keys())]
  operator_index_atom_index = zip(initial_symmetries_operator_indices, initial_symmetries_atom_indices)
  configurations_dict = create_configuration_dict(configurations)
  symmetries = {}
  for operator_index, alt_atom_index in operator_index_atom_index:
      if alt_atom_index == configuration.atom_index:
        continue
      for direction in DIRECTIONS.keys():
        for sign in SIGNS:
          alt_configuration = configurations_dict[alt_atom_index][direction][sign]
          operator = operators[operator_index]
          rotated_configuration = rotate_configuration(configuration, operator) # Obviously don't need to repeat this keep somewhere
          # check equivalence still holds
          rotated_position = rotated_configuration.atoms.get_positions()[configuration.atom_index]
          alt_position = alt_configuration.atoms.get_positions()[alt_configuration.atom_index]
          distance = np.linalg.norm(rotated_position - alt_position)

          if distance <= tol:
            if operator_index in symmetries:    
              if direction in symmetries[operator_index]:
                symmetries[operator_index][direction][sign] = alt_configuration.atom_index
              else:
                symmetries[operator_index][direction] = {}
                symmetries[operator_index][direction][sign] = alt_configuration.atom_index
            else:
              symmetries[operator_index] = {}
              symmetries[operator_index][direction] = {}
              symmetries[operator_index][direction][sign] = alt_configuration.atom_index
  return symmetries

def find_symetries(configurations:List[Configuration], operators, mapping, tol = 1e-2):
  symmetries = {}
  for i, configuration in enumerate(configurations):
    if configuration.atom_index in mapping:
      initial_symmetries = mapping[configuration.atom_index]
      preserved_symmetries = find_symmetries_for_configuration(configuration, configurations, initial_symmetries, operators, tol)
      sign = "+" if configuration.delta > 0 else "-"
      if configuration.atom_index in symmetries:
        if configuration.direction in symmetries[configuration.atom_index]:
          symmetries[configuration.atom_index][configuration.direction][sign] = preserved_symmetries

        else:
          symmetries[configuration.atom_index][configuration.direction] = {}
          symmetries[configuration.atom_index][configuration.direction][sign] = preserved_symmetries
      else:
        symmetries[configuration.atom_index] = {}
        symmetries[configuration.atom_index][configuration.direction] = {}
        symmetries[configuration.atom_index][configuration.direction][sign] = preserved_symmetries
  write_to_json(symmetries, './data/mapping/displaced_atom_mappings.json')
  return symmetries
  