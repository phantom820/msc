from ase.build import bulk
from gpaw import GPAW
import numpy as np

# Start from the primitive 2-atom cell
atoms = bulk('Si', 'diamond', cubic=False)

atoms = bulk('Si', 'diamond', cubic=False)

calc = GPAW(mode='lcao',
            basis='dzp',
            xc='PBE',
            h=0.2,
            symmetry={'point_group': True},
            convergence={"density":1e-5},
            kpts=(3,3, 3),
            txt='./data/si_scf.txt')

atoms.calc = calc

# Compute total energy
E = atoms.get_potential_energy()

# Save symmetry operations
syms = calc.symmetry
np.save('./data/meta/si_symmetries.npy', syms.op_scc)
print(f'Number of symmetry operations: {len(syms.op_scc)}')
print(f"Number of atoms {len(atoms)}")

# Save ground-state wavefunctions (optional)
calc.write('./data/meta/si.gpw')
