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

import netCDF4 as nc
import pymem3dg as dg
import numpy as np


def sizeOf(trajnc: str):
    """Get the number of frame of the netcdf trajectory file
    Args:
        trajNc (str): netcdf trajectory file

    Returns:
        int: number of frames
    """
    with nc.Dataset(trajnc) as ds:
        return ds.groups["Trajectory"].dimensions["frame"].size


def constructSystemByMatrix(trajnc: str, frame: int, parameters: dg.Parameters):
    """Construct pymem3dg System by reading mesh topology and coordinate from netcdf trajectory
    Args:
        trajNc (str): netcdf trajectory file
        frame: frame of the mesh
        parameters: pymem3dg Parameters struct

    Returns:
        int: pymem3dg System object
    """
    with nc.Dataset(trajnc) as ds:
        time = np.array(ds.groups["Trajectory"].variables["time"][frame])
        coordinates = np.array(ds.groups["Trajectory"].variables["coordinates"][frame])
        topology = np.array(ds.groups["Trajectory"].variables["topology"][frame])
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
        return dg.System(topology, coordinates, parameters)


def readMeshByNc(trajNc: str, frame: int):
    """Read topology and coordinate of the mesh from netcdf trajectory

    Args:
        trajNc (str): netcdf trajectory file
        frame (int): frame of the mesh

    Returns:
        ndarray: 2D array of topology matrix
        ndarray: 2D array of vertex coordinate matrix
    """
    with nc.Dataset(trajNc) as ds:
        coordinates = np.array(ds.groups["Trajectory"].variables["coordinates"][frame])
        topology = np.array(ds.groups["Trajectory"].variables["topology"][frame])
        # -1 to dynamically allocate size
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
    return topology, coordinates


def readMeshDataByNc(trajNc: str, frame: int, group: str, variable: str, num_col: int):
    """Read mesh data from netcdf trajectory, example: proteinDensity = readMeshDataByPly(meshList[i], "Trajectory", "proteindensity")

    Args:
        trajNc (str): netcdf trajectory file
        frame (int): frame of the mesh data
        group (str): nc group of the mesh data
        variable (str): mesh data name
        num_col (int): number of columns for reshaping purpose

    Returns:
        ndarray: mesh data in matrix form
    """
    with nc.Dataset(trajNc) as ds:
        data = np.array(ds.groups[group].variables[variable][frame])
        # -1 to dynamically allocate size
        data = np.reshape(data, (-1, num_col))
    return np.squeeze(data)


def readMeshByPly(ply: str):
    """Read topology and coordinate matrix from .ply file

    Args:
        ply (str): .ply file name

    Returns:
        ndarray: 2D array of topology matrix
        ndarray: 2D array of coordinate matrix
    """
    face, vertex = dg.readMesh(ply)
    return face, vertex


def readMeshDataByPly(ply: str, dataType: str, dataName: str):
    """Read mesh data from .ply file, example: proteinDensity = readMeshDataByPly(meshList[i], "vertex", "protein_density")

    Args:
        ply (str): .ply file name
        dataType (str): data type of the mesh data
        dataName (str): name of the mesh data

    Returns:
        ndarray: mesh data
    """
    data = dg.readData(ply, dataType, dataName)
    return data


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
