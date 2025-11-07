import numpy as np


def calculate_distance(ref_position_scaled, position_scaled, epsilon):
  x = ref_position_scaled.copy()
  y = position_scaled.copy()
  x[(1 - x) < epsilon] = 0
  y[(1 - y) < epsilon] = 0
  distance = np.linalg.norm(y - x, axis=1)
  return distance
  