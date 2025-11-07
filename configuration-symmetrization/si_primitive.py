from ase.build import bulk
from gpaw import GPAW
import numpy as np

# Build the primitive diamond structure (2 atoms)
atoms = bulk('Si', 'diamond', cubic=False)

# Attach a GPAW calculator
calc = GPAW(mode='lcao',
            basis='dzp',
            xc='PBE',
            h=0.2,
            symmetry={'point_group': True},
            convergence={'density': 1e-5},
            kpts=(6, 6, 6),
            txt='./data/meta/si_primitive.txt')

atoms.calc = calc

# Compute total energy
E = atoms.get_potential_energy()

# Save symmetry operations
syms = calc.symmetry
np.save('./data/meta/si_primitive_symmetries.npy', syms.op_scc)
print(f'Number of symmetry operations: {len(syms.op_scc)}')

# Save ground-state wavefunctions (optional)
calc.write('./data/meta/si_primitive.gpw')
