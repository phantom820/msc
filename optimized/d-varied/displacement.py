import sys
from ase.build import bulk
from gpaw import GPAW
from gpaw.elph import DisplacementRunner
import numpy as np

# if len(sys.argv) < 2:
#     raise ValueError("Usage: python3 si_displacement.py <h_value>")

h_value = 0.2
d_value = 1e-2

atoms = bulk('Si', 'diamond', cubic=False)

calc = GPAW(
    mode='lcao',
    basis='dzp',
    xc='PBE',
    h=h_value,
    symmetry={'point_group': False},
    convergence={"density": d_value},
    kpts=(3, 3, 3),
    txt=f'./data/si_{d_value}.txt'
)

atoms.calc = calc
atoms.get_potential_energy()

elph = DisplacementRunner(
    atoms=atoms,
    calc=calc,
    supercell=(3, 3, 3),   # << UPDATED
    name='elph',
    calculate_forces=True
)

elph.run()

calc.write(f'./data/si_{d_value}_groundstate.gpw', mode='all')
