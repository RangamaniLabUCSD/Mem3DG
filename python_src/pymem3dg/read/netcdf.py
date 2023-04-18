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
import numpy as np
import numpy.typing as npt
import xarray as xr


def sizeOf(trajnc: str) -> int:
    """Get the number of frame of the netcdf trajectory file
    Args:
        trajNc (str): netcdf trajectory file

    Returns:
        int: number of frames
    """
    with nc.Dataset(trajnc) as ds:
        return ds.groups["Trajectory"].dimensions["frame"].size


def cropNetcdf(trajnc: str) -> None:
    ds_disk = xr.open_dataset(trajnc, group="Trajectory")
    ds_disk = ds_disk.isel(frame=slice(10, 12))
    ds_disk["notablevertex"]
    ds_disk.drop("notablevertex")


def getFaceAndVertexMatrix(
    trajNc: str, frame: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Read topology and coordinate of the mesh from netcdf trajectory

    Args:
        trajNc (str): netcdf trajectory file
        frame (int): frame of the mesh

    Returns:
        tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: 2D array of topology matrix and vertex coordinate matrix
    """
    with nc.Dataset(trajNc) as ds:
        coordinates = np.array(ds.groups["Trajectory"].variables["coordinates"][frame])
        topology = np.array(ds.groups["Trajectory"].variables["topology"][frame])
        # -1 to dynamically allocate size
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
    return topology, coordinates


def getData(
    trajNc: str, frame: int, group: str, variable: str, num_col: int
) -> npt.NDArray[np.float64]:
    """Read mesh data from netcdf trajectory, example: proteinDensity = readMeshDataByPly(meshList[i], "Trajectory", "proteindensity")

    Args:
        trajNc (str): netcdf trajectory file
        frame (int): frame of the mesh data
        group (str): nc group of the mesh data
        variable (str): mesh data name
        num_col (int): number of columns for reshaping purpose

    Returns:
        npt.NDArray[np.float64]: mesh data
    """
    with nc.Dataset(trajNc) as ds:
        data = np.array(ds.groups[group].variables[variable][frame])
        # -1 to dynamically allocate size
        data = np.reshape(data, (-1, num_col))
    return np.squeeze(data)
