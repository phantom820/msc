import json
import numpy as np

def load_json_file(path):
  with open(path, 'r') as f:
    data = json.load(f)
  return data

def potential_from_cache(cache):
  shape = cache["Vt_sG"]["__ndarray__"][0]
  data = np.array(cache["Vt_sG"]["__ndarray__"][2])
  return data.reshape(shape)

def grid_to_frac(shape):
  nx, ny, nz = shape
  grid_coords = np.mgrid[0:nx, 0:ny, 0:nz].reshape(3, -1).T  # shape (nx*ny*nz, 3)
  # frac_coords = grid_coords / np.array([nx, ny, nz])
  frac_coords = grid_coords
  return frac_coords

def verify_potential_values(potential, mapping, shape, threshold = 1e-3):
  flattened_potentital = potential.ravel()
  sym = np.array(mapping["Sym"]).reshape(3,3)
  alt_idx = mapping["M"]
  print("Symetry Matrix")
  print(sym)
  print()

  for i,j in enumerate(mapping["M"]):
    v_i = flattened_potentital[i]
    v_j = flattened_potentital[j]
    print(f"i : {i} V(i) : {v_i} ->  j : {j} V(j) : {v_j} abs diff <= {str(threshold)}: {abs(v_i - v_j) <= threshold}")
    if i == 10:
      break

  diff = np.abs(flattened_potentital - flattened_potentital[alt_idx])
  diff_below_threshold = np.count_nonzero(diff < threshold)
  diff_above_threshold = flattened_potentital.shape[0] - diff_below_threshold
  # print(diff.shape)
  # faulty_idx = np.where(diff > threshold)[0]
  # orig_faulty_idx = np.array(np.unravel_index(faulty_idx, shape))
  # # orig_faulty_idx_frac = grid_to_frac(orig_faulty_idx)
  # print(orig_faulty_idx)
  # print("SS")

  print()
  print(f"Abs diff below threshold count: {diff_below_threshold}")
  print(f"Abs diff above threshold count: {diff_above_threshold}")
  print(f"Maximum abs diff: {diff.max()}")

mappings = load_json_file('./data/mappings.json')
eq_cache = load_json_file('./elph/cache.eq.json')
eq_potential = potential_from_cache(eq_cache)
frac = np.load("./data/frac.npy")
shape = (68,68,68)


# mapping = mappings["2"]
# R = np.array(mapping["Sym"]).reshape(3,3)
# print(R)

# r = np.array([0,0,1])
# r_frac = r/96
# r_frac_p = r_frac@R
# r_grid_p =r_frac_p * shape

# print(r_frac)
# print(r_grid_p)
# a = eq_potential[0, 0, 0, 1]
# b = eq_potential[0, 0, 0, -1]
# print(frac[-1])
# print(eq_potential[0, 0,0, :])

for i, key in enumerate(mappings.keys()):
  verify_potential_values(eq_potential, mappings[key], shape)
  print()
  print("=" * 60)
