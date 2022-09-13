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

import pymem3dg as dg
import numpy as np
import numpy.typing as npt
import pymem3dg.util as dg_util
import copy


def setDefaultMeshProcessor(system: dg.System, lengthScale: float):
    """a commonly used mesh processor

    Args:
        system (dg.System): system handle
        lengthScale (float): lengthScale of the system
    """
    system.meshProcessor.meshMutator.isShiftVertex = True
    system.meshProcessor.meshMutator.flipNonDelaunay = True
    system.meshProcessor.meshMutator.splitLarge = True
    system.meshProcessor.meshMutator.splitFat = True
    system.meshProcessor.meshMutator.splitSkinnyDelaunay = True
    system.meshProcessor.meshMutator.splitCurved = True
    system.meshProcessor.meshMutator.minimumEdgeLength = 0.02 * lengthScale
    system.meshProcessor.meshMutator.maximumEdgeLength = 0.2 * lengthScale
    system.meshProcessor.meshMutator.curvTol = 0.1 / lengthScale
    system.meshProcessor.meshMutator.collapseSkinny = True
    system.meshProcessor.meshMutator.collapseFlat = True
    system.meshProcessor.meshMutator.collapseSmall = True
    system.meshProcessor.meshMutator.targetFaceArea = 0.0003 * lengthScale**2
    system.meshProcessor.meshMutator.isSmoothenMesh = True


def constantOsmoticPressureModel(volume, pressure):
    return pressure, -pressure * volume


def preferredVolumeOsmoticPressureModel(volume, strength, target_volume):
    volume_difference = volume - target_volume
    pressure = -(strength * volume_difference / target_volume / target_volume)
    return pressure, -pressure * volume_difference / 2


def ambientSolutionOsmoticPressureModel(
    volume, enclosed_solute, ambient_concentration, temperature
):
    kBoltzmann = 1.380649e-8 # Boltzmann constant (nanonewton * micrometer / Kelvin)
    i = 1.0 # van't Hoff index
    N = 6.02214076e5 # Avogadro constant (/atto-mol)
    R = kBoltzmann * N # ideal gas constant (nanonewton * micrometer / Kelvin / atto-mol)
    pressure = i * R * temperature * (enclosed_solute / volume - ambient_concentration)
    
    ratio = ambient_concentration * volume / enclosed_solute;
    energy = i * R * temperature * enclosed_solute * (ratio - np.log(ratio) - 1);
    return pressure, energy

def prescribeGeodesicPoteinDensityDistribution(
    time: float,
    vertexMeanCuravtures: npt.NDArray[np.float64],
    geodesicDistance: npt.NDArray[np.float64],
    sharpness: float,
    radius: float,
):
    """form function that prescribe geodesic protein density profile on mesh

    Args:
        time (float): time of the system
        vertexMeanCuravtures (npt.NDArray[np.float64]): vertex mean curvature
        geodesicDistance (npt.NDArray[np.float64]): vertex geodesic curvature
        sharpness (float): sharpness of transition used in tanh function
        radius (float): center of transition of the tanh function

    Returns:
        _type_: _description_
    """
    return dg_util.tanhDistribution(
        x=geodesicDistance, sharpness=sharpness, center=radius
    )


def prescribeGaussianPointForce(
    vertexPositions: npt.NDArray[np.float64],
    vertexDualAreas: list,
    time: float,
    geodesicDistances: list,
    Kf: float,
    std: float,
):
    """form function that generate external force on a single vertex following a Gaussian distribution of geodesic distance

    Args:
        vertexPositions (npt.NDArray[np.float64]): vertex position matrix of the mesh
        vertexDualAreas (list): vertex dual area. Unused
        time (float): time of the simulation. Unused
        geodesicDistances (list): geodesic distance centered at the vertex
        Kf (float): force magnitude coefficient
        std (float): standard deviation of the Gaussian distribution

    Returns:
        (npt.NDArray[np.float64]): vertex force matrix
    """
    direction = dg_util.rowwiseNormalize(vertexPositions)
    magnitude = Kf * dg_util.gaussianDistribution(geodesicDistances, 0, std)
    return dg_util.rowwiseScaling(magnitude, direction)


def prescribePeriodicForceOnCylinder(
    vertexPositions: npt.NDArray[np.float64],
    vertexDualAreas: list,
    time: float,
    geodesicDistance: list,
    Kf: float,
    freq: float,
):
    """prescribe a periodic external force profile on a tubular structure

    Args:
        vertexPositions (npt.NDArray[np.float64]): vertex position matrix of the mesh
        vertexDualAreas (list): vertex dual area. Unused
        time (float): time of the simulation. Unused
        geodesicDistance (list): geodesic distance centered at the vertex
        Kf (float): force magnitude coefficient. Kf = 1 means that the average force = 1, with high = 2 and low = 0
        freq (float): normalized wavenumber of the periodicity. There are #freq cycles along the entire height of the structure

    Returns:
        npt.NDArray[np.float64]: vertex force
    """
    height = np.max(vertexPositions[:, 2]) - np.min(vertexPositions[:, 2])

    # radial unit vector
    direction = copy.deepcopy(vertexPositions)
    direction[:, 2] = 0
    direction = dg_util.rowwiseNormalize(direction)

    magnitude = Kf * (1 + np.sin(freq * 2 * np.pi / height * vertexPositions[:, 2]))

    return dg_util.rowwiseScaling(magnitude, direction)
