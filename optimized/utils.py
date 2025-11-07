import json
import numpy as np
import re
import glob

def write_to_json(data, output_path):
  with open(output_path, 'w') as f:
    json.dump(data, f)

def load_json(path):
  with open(path, 'r') as f:
    data = json.load(f)
  return data

def grid_to_scaled_positions(shape):
  nx, ny, nz = shape
  grid_positions = np.mgrid[0:nx, 0:ny, 0:nz].reshape(3, -1).T  # shape (nx*ny*nz, 3)
  scaled_positions = grid_positions / np.array([nx, ny, nz])
  return scaled_positions

def scaled_positions_to_grid_positions(scaled_position, shape):
  return (scaled_position * np.array(shape)).round().astype(int) % shape

def displacement_key_from_cache_path(cache_path):
  if 'eq.json' in cache_path:
    return 'eq'
  else:
    pattern = re.compile(r'.*cache\.(\d+)([xyz])([+-])\.json')
    m = pattern.match(cache_path)
    if m:
        number, axis, direction = m.groups()
        result = f"{number}|{axis}|{direction}"
        return result
    
def potential_from_cache(cache_path):
  cache = load_json(cache_path)
  shape = cache["Vt_sG"]["__ndarray__"][0]
  data = np.array(cache["Vt_sG"]["__ndarray__"][2])
  return { "displacement_key":  displacement_key_from_cache_path(cache_path), "data" : data.reshape(shape)}

def load_h_varied_potentials(path):
  folders = glob.glob(path)
  results = {}
  for folder in folders:
    h = float(folder.split('/')[-1].split('-')[-1])/100
    results[h] = {'data' : {}}
    cache_paths = glob.glob(f'{folder}/*.json')
    for cache_path in cache_paths:
      if 'info' in cache_path:
        cache_info = load_json(cache_path)
        results[h]['cache_info'] = cache_info
      else:
        potential = potential_from_cache(cache_path)
        results[h]['data'][potential['displacement_key']] = potential
  return results

def load_d_varied_potentials(path):
  folders = glob.glob(path)
  results = {}
  for folder in folders:
    d = float(folder.split('/')[-1].split('-d-')[-1])
    results[d] = {'data' : {}}
    cache_paths = glob.glob(f'{folder}/*.json')
    for cache_path in cache_paths:
      if 'info' in cache_path:
        cache_info = load_json(cache_path)
        results[d]['cache_info'] = cache_info
      elif 'eq' in cache_path:
        potential = potential_from_cache(cache_path)
        results[d]['data'][potential['displacement_key']] = potential
  return results


  


