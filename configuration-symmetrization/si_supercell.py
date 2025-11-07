from ase.build import bulk
from gpaw import GPAW
import numpy as np

# Start from the primitive 2-atom cell
atoms = bulk('Si', 'diamond', cubic=False)

# Make a 3×3×3 supercell (54 atoms)
atoms = atoms.repeat((3, 3, 3))
print(f'Supercell atoms: {len(atoms)}')

# Attach GPAW calculator
# Reduce k-points by the same factor (≈ 6 / 3 = 2)
calc = GPAW(mode='lcao',
            basis='dzp',
            xc='PBE',
            h=0.2,
            symmetry={'point_group': True},
            convergence={'density': 1e-5},
            kpts=(2, 2, 2),
            txt='./data/si_supercell.txt')

atoms.calc = calc

# Compute total energy
E = atoms.get_potential_energy()
print(f'3×3×3 Supercell energy: {E:.6f} eV')

# Save symmetry operations
syms = calc.symmetry
np.save('./data/meta/si_supercell_symmetries.npy', syms.op_scc)
print(f'Number of symmetry operations: {len(syms.op_scc)}')

# Save ground-state wavefunctions (optional)
calc.write('./data/meta/si_supercell.gpw')
