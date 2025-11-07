import numpy as np


def calculate_distance(ref_position_scaled, position_scaled, epsilon = 1e-3):
  x = ref_position_scaled.copy()
  y = position_scaled.copy()
  x[(1 - x) < epsilon] = 0
  y[(1 - y) < epsilon] = 0
  diff = y - x
  if diff.ndim > 1:
    distance = np.linalg.norm(diff, axis=1)
  else:
    distance = np.linalg.norm(diff)
  return distance
  