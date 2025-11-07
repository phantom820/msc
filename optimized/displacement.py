import numpy as np
from typing import List

AXES = {
    "x": 0,
    "y": 1,
    "z": 2
}

AXES_INV = {v: k for k, v in AXES.items()}

class Displacement:
    """Represents an atomic displacement along a Cartesian axis."""
    def __init__(self, axis: str, delta: float):
        if axis not in AXES:
            raise ValueError(f"Invalid axis '{axis}', must be one of {list(AXES)}")
        self.axis = axis
        self.delta = float(delta)
        self.sign = "-" if self.delta < 0 else "+"

    def __repr__(self):
        sign = "+" if self.delta >= 0 else "-"
        return f"<Displacement {self.axis}{sign} | Δ={self.delta}>"

    def key(self, atom_index: int) -> str:
        """Return a compact composite key like '0x+'."""
        sign = "+" if self.delta >= 0 else "-"
        return f"{atom_index}{self.axis}{sign}"
    
    def __str__(self):
        return self.key
    
    def copy(self):
        return Displacement(self.axis, self.delta)

    def __eq__(self, other):
        return isinstance(other, Displacement) and self.axis == other.axis and self.delta == other.delta

def create_displacements(absolute_delta) -> List[Displacement]:
  sgns = [-1, 1]
  displacements = []
  for sgn in sgns:
    for axis in AXES:
      displacements.append(Displacement(axis, sgn * absolute_delta))
  return displacements

def displace_atom(atom_index, atoms, displacement):
  atoms_displaced = atoms.copy()
  cell = atoms_displaced.cell 
  scaled_pos = atoms_displaced.get_scaled_positions()
  delta_vec = np.zeros(3)
  delta_vec[AXES[displacement.axis]] =  displacement.delta  # displacement in Å
  scaled_delta_vec = np.linalg.solve(cell.T, delta_vec)  # fractional displacement
  scaled_pos[atom_index] += scaled_delta_vec
  atoms_displaced.set_scaled_positions(scaled_pos % 1)
  return atoms_displaced