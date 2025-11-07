import numpy as np
from symmetry import symmetrize_equilibrium_potential, symmetrize_displaced_potential
from mapping import create_equilibrium_potential_mapping_table
from scipy.interpolate import RegularGridInterpolator


def cubic_spline_resample(ref_data: np.ndarray, target_shape: tuple[int, int, int]) -> np.ndarray:
  """
  Interpolate a 3D potential grid defined on fractional coordinates [0,1]^3
  to a new grid size using cubic spline interpolation.
  """
  nx, ny, nz = ref_data.shape
  x_ref = np.linspace(0, 1, nx)
  y_ref = np.linspace(0, 1, ny)
  z_ref = np.linspace(0, 1, nz)
  
  interpolator = RegularGridInterpolator(
      (x_ref, y_ref, z_ref),
      ref_data,
      method='cubic',
      bounds_error=False,
      fill_value=None,
  )

  x_tgt = np.linspace(0, 1, target_shape[0])
  y_tgt = np.linspace(0, 1, target_shape[1])
  z_tgt = np.linspace(0, 1, target_shape[2])
  
  # Create 3D mesh of fractional coordinates for the target grid
  X, Y, Z = np.meshgrid(x_tgt, y_tgt, z_tgt, indexing='ij')
  pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)

  resampled = interpolator(pts).reshape(target_shape)
  return resampled

def analyze_potential(potential,symmetrized_potential):
  v = potential[0].ravel()
  v_symmetrized = symmetrized_potential[0].ravel()
  mean_abs_err = np.mean(np.abs(v - v_symmetrized))
  max_abs_diff = np.max(np.abs(v - v_symmetrized))
  mean_rel_err = np.mean(np.abs(v - v_symmetrized) / (np.abs(v) + 1e-12))

  return {
    "mean_abs_err": mean_abs_err,
    "max_abs_diff": max_abs_diff,
    "mean_rel_err": mean_rel_err
  }


def create_d_varied_dataset(d_varied_potentials, eq_potential_mapping):
  d_ref = sorted(d_varied_potentials.keys())[0]
  v_ref = d_varied_potentials[d_ref]['data']['eq']
  data = {'d': [], 'displacement_key': [], 'mean_abs_err' : [], 'max_abs_diff' : [], 'mean_rel_err': []}

  for d in d_varied_potentials:
    if d == d_ref:
      continue

    potentials = d_varied_potentials[d]['data']
    for key in potentials:
      v_sym = symmetrize_equilibrium_potential(potentials[key], eq_potential_mapping)
      analysis = analyze_potential(v_ref['data'], v_sym['data'])
      data['d'].append(d)
      data['displacement_key'].append(key)
      data['mean_abs_err'].append(analysis['mean_abs_err'])
      data['mean_rel_err'].append(analysis['mean_rel_err'])
      data['max_abs_diff'].append(analysis['max_abs_diff'])
  return data


def create_h_varied_dataset(h_varied_potentials, displaced_atom_mapping, eq_potential_mappings):
  data = {'h': [], 'displacement_key': [], 'mean_abs_err' : [], 'max_abs_diff' : [], "mean_rel_err" : []}
  for h in h_varied_potentials:
    potentials = h_varied_potentials[h]['data']
    eq_potential_mapping = eq_potential_mappings[h]
    for key in potentials:
      v = potentials[key]
      if key == 'eq':
        v_sym = symmetrize_equilibrium_potential(potentials[key], eq_potential_mapping)
      else:
        v_sym = symmetrize_displaced_potential(potentials[key], eq_potential_mapping, displaced_atom_mapping)

      analysis = analyze_potential(v['data'], v_sym['data'])
      data['h'].append(h)
      data['displacement_key'].append(key)
      data['mean_abs_err'].append(analysis['mean_abs_err'])
      data['mean_rel_err'].append(analysis['mean_rel_err'])
      data['max_abs_diff'].append(analysis['max_abs_diff'])

  return data


def create_h_varied_dataset_interpolated(h_varied_potentials, displaced_atom_mapping, eq_potential_mappings):
  h_ref = sorted(h_varied_potentials.keys())[0]
  potentials_ref = h_varied_potentials[h_ref]['data']
  data = {'h': [], 'displacement_key': [], 'mean_abs_err' : [], 'max_abs_diff' : [], "mean_rel_err" : []}
  for h in h_varied_potentials:
    if h == h_ref:
      continue
    potentials = h_varied_potentials[h]['data']
    eq_potential_mapping = eq_potential_mappings[h]
    for key in potentials:
   
      if key == 'eq':
        v_sym = symmetrize_equilibrium_potential(potentials[key], eq_potential_mapping)
      else:
        v_sym = symmetrize_displaced_potential(potentials[key], eq_potential_mapping, displaced_atom_mapping)

      v_star = cubic_spline_resample(potentials_ref[key]['data'][0], v_sym['data'][0].shape)
      v = {"displacement_key" : key, "data" : np.array([v_star])}
      analysis = analyze_potential(v['data'], v_sym['data'])
      data['h'].append(h)
      data['displacement_key'].append(key)
      data['mean_abs_err'].append(analysis['mean_abs_err'])
      data['max_abs_diff'].append(analysis['max_abs_diff'])
      data['mean_rel_err'].append(analysis['mean_rel_err'])

  return data