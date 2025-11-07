import numpy as np
from gpaw import GPAW
from ase.build import bulk
from gpaw.symmetry import Symmetry
from ase.geometry import wrap_positions

def displace_atoms(atoms, delta):
	directions = {'x': 0, 'y': 1, 'z': 2}
	signs = {'+': 1, '-': -1}
	cell = atoms.get_cell()
	displacements = {}
	for si_atom in [0, 1]:  # indices of Si atoms in primitive cell
		for dir_name, d in directions.items():
			for sign_name, sgn in signs.items():
					atoms_disp = atoms.copy()
					cell = atoms_disp.cell  # (3,3) matrix
					delta_vec = np.zeros(3)
					delta_vec[d] = sgn * delta  # displacement in Å
					delta_scaled = np.linalg.solve(cell.T, delta_vec)  # fractional displacement

					pos_scaled = atoms_disp.get_scaled_positions()
					pos_scaled[si_atom] += delta_scaled
					# pos_scaled = pos_scaled % 1 # If comment this out no jumps
					atoms_disp.set_scaled_positions(pos_scaled)

					# store in dictionary
					key = (si_atom, dir_name, sign_name)
					displacements[key] = pos_scaled.copy()
	return displacements

# --- Step 2: compute symmetry operations ---
def get_symmetry(atoms):
	id_a = list(atoms.get_atomic_numbers())
	cell_cv = np.array(atoms.get_cell())
	pbc_c = np.array(atoms.get_pbc())
	spos_ac = atoms.get_scaled_positions()

	sym = Symmetry(id_a, cell_cv, pbc_c=pbc_c, tolerance=1e-5)
	sym.analyze(spos_ac)

	return sym

def print_symmetry_table(keys, lookup_matrix, delta):
	n_disp, n_sym = lookup_matrix.shape

	# Print header
	print('Symmetry table')
	print('delta: ', delta)
	headers = ["disp"] + [f"  n={n}  " for n in range(n_sym)]
	for h in headers:
			print(f"{h:<12}", end="")
	print()

	for i in range(len(lookup_matrix)):
		# print(lookup_matrix[i])
		print(f"{i} {keys[i]}  {tuple(lookup_matrix[i][0][:3])} {tuple(lookup_matrix[i][1][:3])} ")
	print()

def print_symmetry_position_pairs(displacements, keys, symmetry_table, cell):

	print('Cell')
	for row in cell:
		print(row)
	print()
	for i in range(len(symmetry_table)):
		symmetries = symmetry_table[i]
		r_atoms_frac = displacements[keys[i]]

		# also get Cartesian
		r_atoms_cart = (r_atoms_frac @ cell)

		for sym in symmetries:
				r_prime_atoms_frac = displacements[sym[:3]]
				r_prime_atoms_cart = (r_prime_atoms_frac @ cell)  # Cartesian

				for j in range(len(r_prime_atoms_frac)):
						print(
								f"{j} {keys[i]} "
								f"frac={r_atoms_frac[j]} cart={r_atoms_cart[j]} "
								f"-> {sym[:3]} "
								f"frac={r_prime_atoms_frac[j]} cart={r_prime_atoms_cart[j]}"
						)
				print()
		print("=" * 200)
		print()

def pretty_print_cartesian(positions):
	print('Original atomic positions')
	header = f"{'Atom':<4}\t{'x':<10}\t{'y':<10}\t{'z':<10}"
	print(header)
	for i, (x, y, z) in enumerate(positions):
		print(f"{i:<4}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}")

def build_symmetry_table_vectorized(displacements, op_scc, ft_sc, threshold = 1e-3):
	keys = list(displacements.keys())
	ndisp = len(keys)
	natoms = displacements[keys[0]].shape[0]

	# 3D array: (num_displacements, num_atoms, 3)
	disp_array = np.array([displacements[k] for k in keys])

	n_sym = len(op_scc)
	lookup_matrix = np.empty((ndisp, n_sym), dtype=object)

	for n, (R, t) in enumerate(zip(op_scc, ft_sc)):
			# Apply symmetry to all displacements at once r' = op_scc x r + ft_c ()
			rotated = (np.einsum('ijk,kl->ijl', disp_array, R.T) + 0 ) # shape (ndisp, natoms, 3)

			# Flatten for distance comparison
			rotated_flat = rotated.reshape(ndisp, natoms*3)
			# rotated[(1 - rotated) <= threshold ] = 0
			disp_flat = disp_array.reshape(ndisp, natoms*3)

			# Distance matrix (rotated vs original)
			dists = np.linalg.norm(rotated_flat[:, None, :] - disp_flat[None, :, :], axis=2)  # (ndisp, ndisp)

			# Find closest match
			closest_idx = np.argmin(dists, axis=1)
			for i, idx in enumerate(closest_idx):
					lookup_matrix[i, n] = keys[idx] + (dists[i, idx],)

	return keys, lookup_matrix

# Load ground state calculation to use same settings
calc = GPAW('Si_groundstate.gpw', parallel={'domain': 1, 'band': 1})
atoms = calc.atoms

delta = 0.1  # displacement in Å
displacements_scaled = displace_atoms(atoms, delta)
positions_scaled = atoms.get_scaled_positions()
symmetries = get_symmetry(atoms)

keys,symmetry_table = build_symmetry_table_vectorized(displacements_scaled, symmetries.op_scc, symmetries.ft_sc)
pretty_print_cartesian(atoms.get_positions())
print()
print_symmetry_table(keys, symmetry_table, delta)
print()
print_symmetry_position_pairs(displacements_scaled, keys, symmetry_table, atoms.get_cell())