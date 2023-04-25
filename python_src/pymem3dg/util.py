# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
#

from scipy.special import sph_harm
import pymem3dg as dg
import numpy as np
import numpy.typing as npt


def rowwiseNormalize(matrix: npt.NDArray[np.float64]):
    """Rowwise normalize the matrix

    Args:
        matrix (ndarray): matrix to be normalized rowwise, dimension = 2

    Returns:
        ndarray: rowwise normalized matrix
    """
    # print("norm:", rowwiseNorm(matrix))
    return rowwiseScaling(rowwiseNorm(matrix) ** (-1), matrix)


def rowwiseNorm(matrix: npt.NDArray[np.float64]):
    """rowwise norm

    Args:
        matrix (ndarray): dimension = 2 or 1

    Returns:
        array: vector of rowwise norm
    """
    if np.shape(matrix)[0] == np.size(matrix):
        return np.abs(matrix)
    else:
        return np.sum(matrix**2, axis=1) ** (0.5)


def rowwiseScaling(scaling: npt.NDArray[np.float64], matrix: npt.NDArray[np.float64]):
    """rowwise scaling

    Args:
        scaling (array): vector of scaling coefficient
        matrix (ndarray): matrix to be scaled, dimension = 2 or 1

    Returns:
        ndarray: scaled matrix
    """
    if np.shape(matrix)[0] == np.size(matrix):
        return matrix * scaling
    else:
        return matrix * scaling[:, None]
        # return matrix * scaling.reshape((-1, 1))


def rowwiseDotProduct(a: npt.NDArray[np.float64], b: npt.NDArray[np.float64]):
    """rowwise dot product

    Args:
        a (ndarray): matrix 1, dimension = 2
        b (ndarray): matrix 2, dimension = 2

    Returns:
        array: vector of rowwise dot product
    """
    return np.sum(a * b, axis=1)


def gaussianDistribution(x: npt.NDArray[np.float64], mean: float, sd: float):
    """Gaussian distribution

    Args:
        x (array): distribution argument
        mean (float): distribution mean
        sd (float): standard deviation

    Returns:
        array: unit gaussian distribution
    """
    prob_density = np.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * (2 * np.pi) ** 0.5)
    return prob_density


def tanhDistribution(x: npt.NDArray[np.float64], sharpness: float, center: float = 0):
    """unit tanh distribution y(x) that is y=1 at x=-inf and y=0 at x=+inf
    Args:
        x (array): distribution argument
        sharpness (float): sharpness of transition
        center (float): location of transition between 1 and 0. Defaults to 0

    Returns:
        array: unit tanh distribution
    """
    return 0.5 * (1.0 + np.tanh(sharpness * (center - x)))


def sphericalHarmonicsPerturbation(
    coordinate: npt.NDArray[np.float64],
    m: int,
    n: int,
    amplitude: float,
    origin: list = [0, 0, 0],
):
    """Apply spherical harmonics to perturb the vertex positions of the mesh in the direction of (position - origin)

    Args:
        coordinate (ndarray): vertex coordinate of the mesh
        m (int): spherical harmonics order
        n (int): spherical harmonics
        amplitude (float): amplitude of perturbation
        origin (_type_, optional): origin of direction vector. Defaults to [0, 0, 0].

    Returns:
        ndarray: perturbed coordinate matrix of the mesh
    """
    x = coordinate[:, 0]
    y = coordinate[:, 1]
    z = coordinate[:, 2]
    unit_vector = rowwiseNormalize(coordinate) - origin
    theta = np.arctan2(y, x) + np.pi
    r = (x**2 + y**2) ** 0.5
    phi = np.arctan(z / (r + 1e-4 * np.max(r))) + np.pi / 2
    harmonics = sph_harm(m, n, theta, phi)
    coordinate = coordinate + amplitude * rowwiseScaling(harmonics.real, unit_vector)
    return coordinate


def smooth(y: list, box_pts: int):
    """Apply moving average to smooth the data

    Args:
        y (array): data to be smoothened
        box_pts (int): box points

    Returns:
        array: smoothened data
    """
    box = np.ones(box_pts) / box_pts
    return np.convolve(y, box, mode="same")

    ## more advanced averaging
    # from scipy.signal import savgol_filter
    # return savgol_filter(y, box_pts+1, 2)


def linearSampling(a_list: list, numPts: int):
    """linear Sampling

    Args:
        a_list (array): list to be sampled
        numPts (int): number of sampling point

    Returns:
        array: sampled list
    """
    indices = np.linspace(0, len(a_list) - 1, num=numPts)
    return np.array([a_list[index] for index in indices])


def femtoJToKBT(energy: float, temperature: float):
    """Unit conversion from femto (1e-15) Joule to kBT

    Args:
        energy (float): energy in Femto Joule
        temperature (float): temperature in Kelvin

    Returns:
        float: energy in kBT
    """
    # temperature = 310
    kBT = 1.38e-23 * temperature * 1e15  # 1e-15 J
    return energy / kBT


def zeroPadding(frame: int, length: int = 6, padding: str = "0"):
    """pad frame index with dummy string in the font

    Args:
        frame (int): frame number
        length (int, optional): number of total string length. Defaults to 6.
        padding (str, optional): padding string. Defaults to "0".

    Returns:
        str: padded string
    """
    out = str(frame)
    while len(out) < length:
        out = padding + out
    return out
