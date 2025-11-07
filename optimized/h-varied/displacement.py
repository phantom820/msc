import sys
from ase.build import bulk
from gpaw import GPAW
from gpaw.elph import DisplacementRunner
import numpy as np

if len(sys.argv) < 2:
    raise ValueError("Usage: python3 si_displacement.py <h_value>")

h_value = float(sys.argv[1])

atoms = bulk('Si', 'diamond', cubic=False)

calc = GPAW(
    mode='lcao',
    basis='dzp',
    xc='PBE',
    h=h_value,
    symmetry={'point_group': False},
    convergence={"density": 1e-4},
    kpts=(3, 3, 3),
    txt=f'./data/si_{h_value}.txt'
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

calc.write(f'./data/si_{h_value}_groundstate.gpw', mode='all')
