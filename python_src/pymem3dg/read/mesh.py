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

import pymem3dg as dg
import numpy.typing as npt
import numpy as np


def getFaceAndVertexMatrix(
    meshName: str,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """read face matrix and vertex matrix from .ply file

    Args:
        meshName (str): mesh file name

    Returns:
        tuple: tuple of face and vertex matrices
    """
    return dg.getFaceAndVertexMatrix(meshName)


def getDataElementName(meshName: str) -> list:
    """Retrieve all richData element name from .ply file. Namely the list of the places where data live in, such as vertex, edge or face

    Args:
        meshName (str): mesh file name

    Returns:
        list: list of places where data lives
    """
    return dg.getDataElementName(meshName)


def getDataPropertyName(meshName: str, elementName: str) -> list:
    """etrive all richData property name from .ply file. Namely the list of the data where data on the particular element, such as vertex, edge or face.

    Args:
        meshName (str): mesh file name
        elementName (str): mesh element, such as vertex, edge or face

    Returns:
        tuple: list of the data name
    """
    return dg.getDataPropertyName(meshName, elementName)


def getData(
    meshName: str, elementName: str, propertyName: str
) -> npt.NDArray[np.float64]:
    """read richData from .ply file

    Args:
        meshName (str): mesh file name
        elementName (str): mesh element, such as vertex, edge or face
        propertyName (str): data name, such as proteinDensity

    Returns:
        npt.NDArray: data as numpy array
    """
    return dg.getData(meshName, elementName, propertyName)


def processSoup(meshName: str) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """process soup data in .ply and return face and vertex matrices.

    Args:
        meshName (str): path to the .ply mesh

    Returns:
        tuple: face and vertex matrices
    """
    return dg.processSoup(meshName)
