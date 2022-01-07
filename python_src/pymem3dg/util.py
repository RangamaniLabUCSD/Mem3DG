import numpy as np
from scipy.special import sph_harm


def rowwise_normalize(matrix):
    """
    Implement a function that normalizes each row of the matrix x (to have unit length).
    
    Argument:
    x -- A numpy matrix of shape (n, m)
    
    Returns:
    x -- The normalized (by row) numpy matrix. You are allowed to modify x.
    """
    
    # Compute x_norm as the norm 2 of x. Use np.linalg.norm(..., ord = 2, axis = ..., keepdims = True)
    matrix_norm = np.linalg.norm(matrix, axis=1, keepdims=False)

    return rowwise_scaling(1/matrix_norm, matrix)

def rowwise_scaling(scaling, matrix):
  if np.shape(matrix)[0] == np.size(matrix):
    return matrix * scaling
  else:
    return matrix * scaling[:, None]

def spherical_harmonics_perturbation(coordinate, m, n, amplitude, origin = None):
  if origin == None:
    origin = [0,0,0]
  x = coordinate[:, 0]
  y = coordinate[:, 1]
  z = coordinate[:, 2]
  unit_vector = rowwise_normalize(coordinate) - origin
  theta = np.arctan2(y, x) + np.pi
  r = (x**2 + y**2)**0.5
  phi = np.arctan(z/r) + np.pi /2 
  harmonics = sph_harm(m, n, theta, phi)
  coordinate = coordinate + amplitude * rowwise_scaling(harmonics.real,unit_vector)
  return coordinate