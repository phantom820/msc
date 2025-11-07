
import numpy as np
import json
import glob

def load_symmetries(path):
  return np.load(path)

def apply_sym(frac, R, shape):
  # check mltiply order
  # new_frac = (R @ frac) % 1.0  # wrap back into [0,1)
  new_frac = (frac @ R) % shape
  return new_frac

def frac_to_grid(frac, shape):
  idx = frac % shape
  return tuple(idx)

def grid_to_frac(shape):
  nx, ny, nz = shape
  grid_coords = np.mgrid[0:nx, 0:ny, 0:nz].reshape(3, -1).T  # shape (nx*ny*nz, 3)
  # frac_coords = grid_coords / np.array([nx, ny, nz])
  frac_coords = grid_coords
  return frac_coords


def valdiate_mapping(mapping):
  keys = mapping.keys()
  unique, counts = np.unique(list(mapping.values()), return_counts = True)
  
  assert len(keys) == unique.shape[0], "Not all indices appear in mapping"
  assert np.all(counts == 1) == True, "Duplicate index exists in mapping"

def build_mapping(frac, R, shape):
  mapping = {}
  for i, r in enumerate(frac):
    r_p = apply_sym(r, R, shape)
    idx_p = frac_to_grid(r_p, shape)
    j = int(np.ravel_multi_index(idx_p, shape))
    mapping[i] = j
  return mapping

def build_all_mappings(frac, symmetries, shape):
  mappings = {}
  for i, R in enumerate(symmetries):
    print(f"Building mapping for symmetry: {i}")
    print(R)
    mapping = build_mapping(frac, R, shape)
    valdiate_mapping(mapping)
    mappings[i] = {
      'Sym' : R.ravel().tolist(),
      'Sym_shape': R.shape,
      # 'M': [{k:v} for k,v in mapping.items()]
      'M': list(mapping.values())  # Store better
    }
    print("=" * 60)
  return mappings

def load_eq_potential_h_metadata(path):
  files = glob.glob(path)
  metadata = []
  for file in files:
    with open(file, "r") as f:
      cache = json.load(f)
      shape = cache["Vt_sG"]["__ndarray__"][0][1:]
      h = int(file.split("/")[-1].split(".")[1].split("-")[1])
      metadata.append({"shape" : shape, "h" : h})
  return metadata

def build_h_mappings():
  h_metadata = load_eq_potential_h_metadata("./data/h/eq-potential/*.json")
  for meta in h_metadata:
    shape = meta["shape"]
    filename = f"./data/h/eq-potential-mapping/mapping.h-{meta["h"]}.json"
    frac = grid_to_frac(shape)
    print("Processing metadata:", str(meta))
    all_mappings = build_all_mappings(frac, symmetries, shape)
    with open(filename, 'w') as f:
      json.dump(all_mappings, f)
    print("*" * 60)
    print()



def load_eq_potential_d_metadata(path):
  files = glob.glob(path)
  metadata = []
  for file in files:
    with open(file, "r") as f:
      cache = json.load(f)
      shape = cache["Vt_sG"]["__ndarray__"][0][1:]
      d = float(file.split("cache.eq.")[1].split(".json")[0])
      metadata.append({"shape" : shape, "d" : d})
  return metadata
  
def build_d_mappings(): 
  d_metadata = load_eq_potential_d_metadata("./data/d/eq-potential/*.json")
  for meta in d_metadata:
    shape = meta["shape"]
    filename = f"./data/d/eq-potential-mapping/mapping.d-{meta["d"]}.json"
    frac = grid_to_frac(shape)
    print("Processing metadata:", str(meta))
    all_mappings = build_all_mappings(frac, symmetries, shape)
    with open(filename, 'w') as f:
      json.dump(all_mappings, f)
    print("*" * 60)
    print()



symmetries = load_symmetries('./data/si_symmetries.npy')
build_d_mappings()
build_h_mappings()