import numpy as np
from gpaw import GPAW
from displacement import create_displacements
from configuration import create_configurations
from symmetry import find_symetries, symmetrize_equilibrium_potential, symmetrize_displaced_potential
from mapping import create_equilibrium_atom_mapping_table, create_equilibrium_potential_mapping_table
from utils import potential_from_cache, write_to_json, load_json, load_h_varied_potentials, load_d_varied_potentials
from analysis import analyze_potential, create_h_varied_dataset, create_h_varied_dataset_interpolated, create_d_varied_dataset
import pandas as pd

name = "si"
calc = GPAW(f'./data/meta/{name}.gpw', parallel={'domain': 1, 'band': 1})
operators = np.load(f'./data/meta/{name}_symmetries.npy')
atoms = calc.atoms


# data = analyze_potential(eq_potential, symmetrized_eq_potential)
# data2 = analyze_potential(displaced_potential["data"], symmetrized_displaced_potential["data"])

h_varied_potentials = load_h_varied_potentials('./h-varied/elph-*')
d_varied_potentials = load_d_varied_potentials('./d-varied/elph-*')

h_ref = sorted(h_varied_potentials.keys())[0]
d_ref = sorted(d_varied_potentials.keys())[0]
potential_cache_info_ref = h_varied_potentials[h_ref]['cache_info']
delta = potential_cache_info_ref["delta"]
eq_potential_ref =  h_varied_potentials[h_ref]['data']['eq']['data']

eq_atom_mapping = create_equilibrium_atom_mapping_table(atoms, operators)

h_eq_potential_mapping = { 
  k : create_equilibrium_potential_mapping_table(operators, h_varied_potentials[k]['data']['eq']['data'].shape[1:]) for k in h_varied_potentials}
d_eq_potential_mapping = create_equilibrium_potential_mapping_table(operators, d_varied_potentials[d_ref]['data']['eq']['data'].shape[1:])

displacements = create_displacements(delta)
configurations = create_configurations(atoms, displacements)
displaced_atom_mapping = find_symetries(configurations, operators, eq_atom_mapping)

d_data = create_d_varied_dataset(d_varied_potentials, d_eq_potential_mapping)
df_d = pd.DataFrame(d_data)
df_d.to_csv('./data/d_varied_potential_symmetrization.csv', index=False)

data = create_h_varied_dataset(h_varied_potentials, displaced_atom_mapping, h_eq_potential_mapping)
df = pd.DataFrame(data)
df.to_csv('./data/h_varied_potential_symmetrization.csv', index=False)

data_interpolated = create_h_varied_dataset_interpolated(h_varied_potentials, displaced_atom_mapping, h_eq_potential_mapping)
df_interopolated = pd.DataFrame(data_interpolated)
df_interopolated.to_csv('./data/h_varied_potential_symmetrization_interpolated.csv', index=False)