import numpy as np

DIRECTIONS = {
  'x' : 0,
  'y' : 1,
  'z' : 2
}

SIGNS = ["+", "-"]

class Displacement:
  def __init__(self, direction, delta):
    self.direction = direction
    self.delta = delta
	
  def __str__(self):
    return f"({self.direction}, {self.delta})"
  
def get_all_displacements(absolute_delta):
	sgns = [-1, 1]
	displacements = []
	for sgn in sgns:
		for direction in DIRECTIONS:
			displacements.append(Displacement(direction, sgn * absolute_delta))
	return displacements

def displace_atom(atom_index, atoms, displacement):
	atoms_disp = atoms.copy()
	cell = atoms_disp.cell 
	scaled_pos = atoms_disp.get_scaled_positions()
	delta_vec = np.zeros(3)
	delta_vec[DIRECTIONS[displacement.direction]] =  displacement.delta  # displacement in Å
	scaled_delta_vec = np.linalg.solve(cell.T, delta_vec)  # fractional displacement
	scaled_pos[atom_index] += scaled_delta_vec
	atoms_disp.set_scaled_positions(scaled_pos % 1)
	return atoms_disp
