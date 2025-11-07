from configuration import Configuration


def rotate_configuration(configuration:Configuration, rotation) -> Configuration:
  rotated_atoms = rotate_atom(configuration.atom_index, configuration.atoms.copy(), rotation)
  return Configuration(configuration.atom_index, rotated_atoms, configuration.displacement)

def rotate_atom(atom_index, atoms, rotation):
  scaled_positions = atoms.get_scaled_positions()
  rotated_scaled_position = (scaled_positions[atom_index] @ rotation) % 1
  rotated_atoms = atoms.copy()
  rotated_scaled_positions = rotated_atoms.get_scaled_positions()
  rotated_scaled_positions[atom_index] = rotated_scaled_position
  rotated_atoms.set_scaled_positions(rotated_scaled_positions)
  return rotated_atoms

def rotate_scaled_positions(scaled_positions, rotation):
  return (scaled_positions @ rotation) % 1