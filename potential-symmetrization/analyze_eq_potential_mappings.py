import numpy as np
import glob
import json
from scipy.interpolate import RegularGridInterpolator


def load_json_file(path):
  with open(path, 'r') as f:
    data = json.load(f)
  return data

def potential_from_cache(path):
  cache = load_json_file(path)
  shape = cache["Vt_sG"]["__ndarray__"][0]
  data = np.array(cache["Vt_sG"]["__ndarray__"][2])
  return data.reshape(shape)

def symetrize_potential(potential, sym):
  flattened_potentital = potential.ravel()
  alt_idx = sym["M"]
  return flattened_potentital[alt_idx]

def symetrize_and_aggregate(potential, mapping):
  symetrized_potential = np.zeros(potential.ravel().shape)
  n = 0
  for i, key in enumerate(mapping.keys()):
    v = symetrize_potential(potential, mapping[key])
    symetrized_potential = symetrized_potential + v
    n = n + 1
  return symetrized_potential/n


def get_aggregate_h_symetrized_potentials(h_data):
  data = []
  for x in h_data:
    h = x["h"]
    potential = x["potential"]
    mapping = x["mapping"]
    sym_potential = symetrize_and_aggregate(potential, mapping)
    data.append({"h" : h, "v" : potential.ravel(), "v_sym" : sym_potential, "shape" : potential.shape})
  return data


def get_aggregate_d_symetrized_potentials(h_data):
  data = []
  for x in h_data:
    d = x["d"]
    potential = x["potential"]
    mapping = x["mapping"]
    sym_potential = symetrize_and_aggregate(potential, mapping)
    data.append({"d" : d, "v" : potential.ravel(), "v_sym" : sym_potential, "shape" : potential.shape})
  return data
  

def load_h_varied_data(path):
  files = sorted(glob.glob(path))
  mappings = []
  for file in files:
    h = int(file.split('/')[-1].split(".")[1].split("-")[1])
    mapping = load_json_file(file)
    potential = potential_from_cache(f"./data/h/eq-potential/cache.h-{h}.json")
    mappings.append({"h":h/100, "mapping": mapping, "potential": potential })
  return mappings



def load_d_varied_data(path):
  files = glob.glob(path)
  mappings = []
  ds = []
  for file in files:
    d = float(file.split("mapping.d-")[1].split(".json")[0])
    ds.append(d)
    mapping = load_json_file(file)
    potential = potential_from_cache(f"./data/d/eq-potential/cache.eq.{d}.json")
    mappings.append({"d":d, "mapping": mapping, "potential": potential })
  idx = np.argsort(ds)
  return [mappings[i] for i in idx]


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


def create_data(data):
  h = []
  mean_abs_err = []
  max_abs_diff = []
  for point in data:
    h.append(point["h"])
    max_abs_diff.append(np.max(np.abs(point["v"] - point["v_sym"])))
    mean_abs_err.append(np.mean(np.abs(point["v"] - point["v_sym"])))

  return {
    "h": h,
    "mean_abs_err": np.array(mean_abs_err),
    "max_abs_diff": np.array(max_abs_diff)
  }

def create_d_data(data):
  d = []
  mean_abs_err = []
  max_abs_diff = []
  for point in data:
    d.append(point["d"])
    max_abs_diff.append(np.max(np.abs(point["v"] - point["v_sym"])))
    mean_abs_err.append(np.mean(np.abs(point["v"] - point["v_sym"])))

  return {
    "d": d,
    "mean_abs_err": np.array(mean_abs_err),
    "max_abs_diff": np.array(max_abs_diff)
  }

def create_data_interpolated(ref_data,  data):
  h = []
  ref_sum = []
  sym_sum = []
  mean_abs_err = []
  max_abs_diff = []
  ref_v_shape = ref_data["shape"]
  ref_v = ref_data["v"].reshape(ref_v_shape)
  for point in data:
    target_shape = point["shape"][1:]
    ref_v_interpolated = cubic_spline_resample(ref_v[0], target_shape).ravel()
    h.append(point["h"])
    max_abs_diff.append(np.max(np.abs(ref_v_interpolated - point["v_sym"])))
    mean_abs_err.append(np.mean(np.abs(ref_v_interpolated - point["v_sym"])))
    ref_sum.append(ref_v_interpolated.sum()/ ref_v_interpolated.shape[0])
    sym_sum.append(point["v_sym"].sum()/ point["v_sym"].shape[0] )

  return {
    "h": h,
    "ref_sum": np.array(ref_sum),
    "sym_sum": np.array(sym_sum),
    "mean_abs_err": np.array(mean_abs_err),
    "max_abs_diff": np.array(max_abs_diff)
  }



if __name__ == "__main__":
    
  # h_varied_mappings = load_h_varied_data("./data/h/eq-potential-mapping/*.json")
  # raw_data = get_aggregate_h_symetrized_potentials(h_varied_mappings)
  # ref_data = raw_data[0]

  # # data = create_data(raw_data)
  # data_interpolated = create_data_interpolated(ref_data, raw_data[1:])
  # print(data_interpolated)

  d_varied_mappings = load_d_varied_data("./data/d/eq-potential-mapping/*.json")
  raw_d_data = get_aggregate_d_symetrized_potentials(d_varied_mappings)
  x = create_d_data(raw_d_data)
  print(x)
  # print(ref_data["shape"])
  # h_varied_mappings = h_varied_mappings[1:]
  # x = get_aggregate_h_symetrized_potentials(h_varied_mappings)

  # for point in x:
  #   # point = x[0]ß
  #   h = point["h"] 
  #   max_diff = np.max(np.abs(point["v"] - point["v_sym"]))
  #   mean_abs = np.mean(np.abs(point["v"] - point["v_sym"]))

  #   print(h)
  #   print(max_diff)
  #   print()