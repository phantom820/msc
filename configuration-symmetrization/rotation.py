import numpy as np
from gpaw import GPAW
from ase.build import bulk
from gpaw.symmetry import Symmetry
from utils import write_to_json
from displacement import get_all_displacements, displace_atom
from configuration import Configuration
from ase.geometry import wrap_positions


def rotate_configuration(configuration, rotation):
  rotated_atoms = rotate_atom(configuration.atom_index, configuration.atoms.copy(), rotation)

  return Configuration(configuration.atom_index, rotated_atoms, configuration.direction, configuration.delta)

def rotate_atom(atom_index, atoms, rotation):
  scaled_positions = atoms.get_scaled_positions()
  # rotated_scaled_position = wrap_positions([scaled_positions[atom_index] @ rotation], atoms.cell)[0]
  rotated_scaled_position = (scaled_positions[atom_index] @ rotation) % 1

  rotated_atoms = atoms.copy()
  rotated_scaled_positions = rotated_atoms.get_scaled_positions()
  rotated_scaled_positions[atom_index] = rotated_scaled_position
  rotated_atoms.set_scaled_positions(rotated_scaled_positions)
  return rotated_atoms

def rotate_atoms(atoms, rotation):
  scaled_positions = atoms.get_scaled_positions()
  rotated_scaled_positions = scaled_positions @ rotation
  rotated_atoms = atoms.copy()
  rotated_atoms.set_scaled_positions(rotated_scaled_positions)
  return rotated_atoms