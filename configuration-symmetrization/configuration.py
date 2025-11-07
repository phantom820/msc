

from displacement import displace_atom

class Configuration:
	def __init__(self, atom_index, atoms, direction, delta):
		self.atom_index = atom_index
		self.atoms = atoms
		self.direction = direction
		self.delta = delta

	def __str__(self):
		return f"({self.atom_index}, {self.direction}, {self.delta})"
	
def create_configurations(atoms, displacements):
  configurations = []
  for i, atom in enumerate(atoms):
    for displacement in displacements:
      displaced_atoms = displace_atom(i, atoms, displacement)
      configuration = Configuration(i, displaced_atoms, displacement.direction, displacement.delta)
      configurations.append(configuration)
  return configurations