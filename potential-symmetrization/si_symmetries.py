from ase.build import bulk
from gpaw import GPAW
import numpy as np

atoms = bulk('Si', 'diamond', cubic=False)

# Not min h -> 0.14
calc = GPAW(mode='lcao',
            basis='dzp',
            xc='PBE',
            h=0.2,
            symmetry={'point_group': True},
            convergence={"density":1e-5},
            kpts=(6,6, 6),
            txt='./data/si_scf.txt')

atoms.calc = calc

atoms.get_potential_energy()

syms = calc.symmetry
np.save('./data/si_symmetries.npy', syms.op_scc)
print(syms)
