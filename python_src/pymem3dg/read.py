import netCDF4 as nc
import pymem3dg as dg
import numpy as np


def sizeOf(trajnc):
    with nc.Dataset(trajnc) as ds:
        return ds.groups["Trajectory"].dimensions['frame'].size


def constructSystemByMatrix(trajnc, frame, parameters):
    with nc.Dataset(trajnc) as ds:
        time = np.array(ds.groups['Trajectory'].variables['time'][frame])
        coordinates = np.array(
            ds.groups['Trajectory'].variables['coordinates'][frame])
        topology = np.array(
            ds.groups['Trajectory'].variables['topology'][frame])
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
        return dg.System(topology, coordinates, parameters)


def readMeshByNc(trajNc, frame):
    with nc.Dataset(trajNc) as ds:
        coordinates = np.array(
            ds.groups['Trajectory'].variables['coordinates'][frame])
        topology = np.array(
            ds.groups['Trajectory'].variables['topology'][frame])
        # -1 to dynamically allocate size
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
    return topology, coordinates


def readMeshDataByNc(trajNc, frame, group, variable, num_col):
    """ example: proteinDensity = readMeshDataByPly(meshList[i], "Trajectory", "proteindensity") """
    with nc.Dataset(trajNc) as ds:
        data = np.array(
            ds.groups[group].variables[variable][frame])
        # -1 to dynamically allocate size
        data = np.reshape(data, (-1, num_col))
    return np.squeeze(data)


def readMeshByPly(ply):
    face, vertex = dg.readMesh(ply)
    return face, vertex


def readMeshDataByPly(ply, dataType, dataName):
    """ example: proteinDensity = readMeshDataByPly(meshList[i], "vertex", "protein_density") """
    data = dg.readData(ply,
                       dataType, dataName)
    return data
