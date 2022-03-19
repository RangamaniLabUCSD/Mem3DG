from scipy.special import sph_harm
import pymem3dg as dg
import numpy as np


def rowwiseNormalize(matrix):
    # print("norm:", rowwiseNorm(matrix))
    return rowwiseScaling(rowwiseNorm(matrix)**(-1), matrix)


def rowwiseNorm(matrix):
    if np.shape(matrix)[0] == np.size(matrix):
        return np.abs(matrix)
    else:
        return np.sum(matrix**2, axis=1)**(0.5)


def rowwiseScaling(scaling, matrix):
    if np.shape(matrix)[0] == np.size(matrix):
        return matrix * scaling
    else:
        return matrix * scaling[:, None]


def rowwiseDotProduct(a, b):
    return np.sum(a*b, axis=1)


def spherical_harmonics_perturbation(coordinate, m, n, amplitude, origin=None):
    if origin == None:
        origin = [0, 0, 0]
    x = coordinate[:, 0]
    y = coordinate[:, 1]
    z = coordinate[:, 2]
    unit_vector = rowwiseNormalize(coordinate) - origin
    theta = np.arctan2(y, x) + np.pi
    r = (x**2 + y**2)**0.5
    phi = np.arctan(z/r) + np.pi / 2
    harmonics = sph_harm(m, n, theta, phi)
    coordinate = coordinate + amplitude * \
        rowwiseScaling(harmonics.real, unit_vector)
    return coordinate

def smooth(y, box_pts):
    """moving average"""
    box = np.ones(box_pts) / box_pts
    return np.convolve(y, box, mode="same")

    """ more advanced """
    # from scipy.signal import savgol_filter
    # return savgol_filter(y, box_pts+1, 2)
