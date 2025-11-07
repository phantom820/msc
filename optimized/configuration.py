from displacement import Displacement, displace_atom
from typing import List

class Configuration:
    """
    Represents a configuration where one atom is displaced and all others remain unchanged.
    """
    def __init__(self, atom_index: int, atoms, displacement: Displacement):
        self.atom_index = atom_index
        self.atoms = atoms  # can be ASE Atoms object or np.array(N,3)
        self.displacement = displacement
        self.key =  f"{self.atom_index}|{self.displacement.axis}|{self.displacement.sign}"
    
    @staticmethod
    def create_key(atom_index, axis, sign):
       return f"{atom_index}|{axis}|{sign}"
    
    def __str__(self):
       return self.key
    
def create_configurations(atoms, displacements:List[Displacement]) -> List[Configuration]:
  configurations = []
  for atom_index in range(len(atoms)):
    for displacement in displacements:
      displaced_atoms = displace_atom(atom_index, atoms, displacement)
      configuration = Configuration(atom_index, displaced_atoms, displacement)
      configurations.append(configuration)
  return configurations