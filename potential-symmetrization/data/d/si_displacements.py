from ase.build import bulk
from gpaw import GPAW
from gpaw.elph import DisplacementRunner
import numpy as np

atoms = bulk('Si', 'diamond', cubic=False)

# Not min h -> 0.14
# h = 0.14 , step size 0.01 till 0.20 h was fixed at 1e-5
# conv array([1.000e-06, 2.575e-05, 5.050e-05, 7.525e-05, 1.000e-04]) h was fixed at 0.2
calc = GPAW(mode='lcao',
            basis='dzp',
            xc='PBE',
            h=0.12,
            symmetry={'point_group': False},
            convergence={"density":1e-5},
            kpts=(6,6, 6),
            txt='./data/si_scf.txt')

atoms.calc = calc

atoms.get_potential_energy()

elph = DisplacementRunner(atoms=atoms,
                          calc=calc,
                          supercell=(3, 3, 3),
                          name='elph',
                          calculate_forces=True)

elph.run()

calc.write('./data/si_groundstate.gpw', mode='all')
# converge 1e-4 to 1e-6
# h 0.14 -> 0.2 practical 