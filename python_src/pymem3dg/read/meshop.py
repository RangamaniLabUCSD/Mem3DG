#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2024- The Mem3DG Authors
# and the project initiators Cuncheng Zhu, Christopher T. Lee, and
# Padmini Rangamani.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Please help us support Mem3DG development by citing the research
# papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
# for more information.

from re import sub
import pymem3dg as dg
import numpy as np
import numpy.typing as npt


def getCylinder(
    radius: float,
    radialSubdivision: int,
    axialSubdivision: int,
    frequency: float = 1,
    amplitude: float = 0,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Generate cylinderal mesh

    Args:
        radius (float): radius of the tube
        radialSubdivision (int): number of subdivision in the radial direction
        axialSubdivision (int): # subdivision in the axial direction
        frequency (float, optional): normalized wavenumber. Defaults to 1.
        amplitude (float, optional): amplitude of sinusoidal perturbation. Defaults to 0.

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrices
    """
    return dg.getCylinder(
        radius, radialSubdivision, axialSubdivision, frequency, amplitude
    )


def getIcosphere(
    radius: float, subdivision: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Generate spherical mesh

    Args:
        radius (float): radius of the sphere
        subdivision (int): # of subdivision

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrices
    """
    return dg.getIcosphere(radius, subdivision)


def getHexagon(
    radius: float, subdivision: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """generate hexagonal mesh

    Args:
        radius (float): radius of the hexagonal patch
        subdivision (int): # subdivision

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrices
    """
    return dg.getHexagon(radius, subdivision)


def getTetrahedron() -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """generate tetrahedron mesh

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrices
    """
    return dg.getTetrahedron()


def getDiamond(
    dihedral: float,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """generate diamond shaped mesh

    Args:
        dihedral (float): dihedral angle of the diamond

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrices
    """
    return dg.getDiamond(dihedral)


def getFaceSurfacePointClosestToEmbeddedCoordinate(
    faceMatrix: npt.NDArray[np.int64],
    vertexMatrix: npt.NDArray[np.float64],
    embeddedCoordinate: list,
    filter: list = None,
    accountedCoordinate: list = [True, True, True],
) -> tuple[int, list]:
    """Find the face index and the barycentric coordinate of the face surface
    point closest to a embedded coordinate in Euclidean distance

       Args:
           faceMatrix (npt.NDArray[np.int64]): face topology matrix of the mesh (F x 3)
           vertexMatrix (npt.NDArray[np.float64]): vertex position matrix of the mesh (V x 3)
           embeddedCoordinate (list): the target embedded coordinate expressed in 3D cartesian coordinate
            filter (list): filter Limit the scope of search to a subset of vertices, Defaults to all true
           accountedCoordinate (list, optional): array to identify the accounted coordinate. For example, it finds closest vertex in the x-y plane when set to be {true, true, false}. Defaults to {true, true, true}

       Returns:
           tuple[int, list]: tuple of face index and the barycentric coordinate
    """
    if filter == None:
        filter = [True] * np.shape(vertexMatrix)[0]
    return dg.getFaceSurfacePointClosestToEmbeddedCoordinate(
        faceMatrix=faceMatrix,
        vertexMatrix=vertexMatrix,
        embeddedCoordinate=embeddedCoordinate,
        filter=filter,
        accountedCoordinate=accountedCoordinate,
    )


def getVertexClosestToEmbeddedCoordinate(
    vertexMatrix: npt.NDArray[np.float64],
    embeddedCoordinate: list,
    filter: list = None,
    accountedCoordinate: list = [True, True, True],
) -> tuple[int, list]:
    """Find the vertex index closest to a embedded coordinate in Euclidean distance

    Args:
        vertexMatrix (npt.NDArray[np.float64]): vertex position matrix of the mesh (V x 3)
        embeddedCoordinate (list): the target embedded coordinate expressed in 3D cartesian coordinate
        filter (list): filter Limit the scope of search to a subset of vertices, Defaults to all true
        accountedCoordinate (list, optional): array to identify the accounted coordinate. For example, it finds closest vertex in the x-y plane when set to be {true, true, false}. Defaults to {true, true, true}

    Returns:
        tuple[int, list]: tuple of face index and the barycentric coordinate
    """
    if filter == None:
        filter = [True] * np.shape(vertexMatrix)[0]
    return dg.getVertexClosestToEmbeddedCoordinate(
        vertexMatrix=vertexMatrix,
        embeddedCoordinate=embeddedCoordinate,
        filter=filter,
        accountedCoordinate=accountedCoordinate,
    )


def getVertexFurthestFromBoundary(
    faceMatrix: npt.NDArray[np.int64], vertexMatrix: npt.NDArray[np.float64]
) -> int:
    """Find the index of the vertex furthest away from the boundaries

    Args:
        faceMatrix (npt.NDArray[np.int64]): face matrix
        vertexMatrix (npt.NDArray[np.float64]): vertex matrix

    Returns:
        int: vertex index
    """
    return dg.getVertexFurthestFromBoundary(faceMatrix, vertexMatrix)


def linearSubdivide(
    face: npt.NDArray[np.int64], vertex: npt.NDArray[np.float64], nSub: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """subdivide the mesh with linear interpolation

    Args:
        face (npt.NDArray[np.int64]): face matrix
        vertex (npt.NDArray[np.float64]): vertex matrix
        nSub (int): # subdivision

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrix
    """
    return dg.linearSubdivide(face=face, vertex=vertex, nSub=nSub)


def loopSubdivide(
    face: npt.NDArray[np.int64], vertex: npt.NDArray[np.float64], nSub: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """subdivide the mesh with loop scheme

    Args:
        face (npt.NDArray[np.int64]): face matrix
        vertex (npt.NDArray[np.float64]): vertex matrix
        nSub (int): # subdivision

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: face and vertex matrix
    """
    return dg.loopSubdivide(face=face, vertex=vertex, nSub=nSub)
