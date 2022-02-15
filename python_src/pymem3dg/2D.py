import autograd.numpy as np
from autograd import elementwise_grad as egrad
from autograd import grad
from autograd import value_and_grad
from ._core import *


class Parameters:
    Kb = 0
    Kbc = 0
    Ksg = 0
    At = 0
    epsilon = 0
    Kv = 0
    Vt = 0

    def __init__(self) -> None:
        pass


class Geometry:
    nVertex = 0
    nEdge = 0
    ifClosed = True

    vertexTurningAngles = np.array([])
    vertexPositions = np.array([])

    edgeLengths = np.array([])
    edgeAbsoluteAngles = np.array([])

    def __init__(self, vertexPositions, ifClosed):
        self.vertexPositions = vertexPositions
        self.populate_geometry(ifClosed)

    def populate_geometry(self, ifClosed):
        self.ifClosed = ifClosed
        x = self.vertexPositions[:, 0]
        y = self.vertexPositions[:, 1]
        self.nVertex = self.vertexPositions.shape[0]
        if ifClosed:
            self.nEdge = self.geometry.nVertex
            dx = np.diff(np.append(x, x[0]))  # n+1-1
            dy = np.diff(np.append(y, y[0]))
            self.edgeLengths = np.sqrt(dx**2 + dy**2)  # n
            self.edgeAbsoluteAngles = np.arctan2(dy, dx)  # n
            self.vertexTurningAngles = np.diff(np.append(self.edgeAbsoluteAngles,
                                                         self.edgeAbsoluteAngles[0]))  # n
        else:
            raise Exception("Not implemented yet!")

    def check_dimensionality(self):
        if True:
            return 0


class System:
    class Force:
        bendingForce = np.array([])
        surfaceForce = np.array([])
        pressureForce = np.array([])
        adsorptionForce = np.array([])

        def __init__(self) -> None:
            pass

    class Energy:
        totalEnergy = 0
        bendingEnergy = 0
        surfaceEnergy = 0
        pressureEnergy = 0
        adsorptionEnergy = 0

        def __init__(self) -> None:
            pass

    class GeometricMeasurements:
        vertexArclengthForms = np.array([])
        edgeAreaForms = np.array([])
        # vertexVolumeForms = np.array([])

        vertezCurvatures = np.array([])

    def __init__(self, geometry, parameters):
        self.energy = self.Energy()
        self.force = self.Force()
        self.geometricMeasurements = self.GeometricMeasurements()

        self.geometry = geometry
        self.parameters = parameters
        self.compute_totalEnergy(
            self.join(self.geometry.edgeLengths, self.geometry.vertexTurningAngles))

    def compute_totalEnergy(self, vertexPosition_flattened):
        x, y = self.split(vertexPosition_flattened)
        self.geometry.nEdge = self.geometry.nVertex
        dx = np.diff(np.append(x, x[0]))  # n+1-1
        dy = np.diff(np.append(y, y[0]))
        self.geometry.edgeLengths = np.sqrt(dx**2 + dy**2)  # n
        self.edgeAbsoluteAngles = np.arctan2(dy, dx)  # n
        self.vertexTurningAngles = np.diff(np.append(self.edgeAbsoluteAngles,
                                                      self.edgeAbsoluteAngles[0]))  # n

        self.geometricMeasurements.vertexArclengthForms = (
            np.append(edgeLengths[-1], edgeLengths)[:self.geometry.nEdge] + edgeLengths) * 0.5

        self.geometricMeasurements.vertexCurvatures = vertexTurningAngles / \
            self.geometricMeasurements.vertexArclengthForms

        self.geometricMeasurements.edgeAreaForms = 2 * np.pi * radius * self.arclength_v

        self.volume_v = np.pi * radius_e**2 * dh

        self.energy.bendingEnergy = self.parameters.Kb * \
            np.sum((self.meanCurvature_v - 10)**2 * self.area_v)

        self.energy.surfaceEnergy = self.parameters.Ksg * \
            (np.sum(self.area_v) - self.parameters.At)**2

        self.energy.pressureEnergy = self.parameters.Kv * \
            (np.sum(self.volume_v) - self.parameters.Vt)**2

        targetLength = ds / ds.shape[0]
        regularization = np.linalg.norm(ds - targetLength)**2
        return self.energy.bendingEnergy + self.energy.surfaceEnergy + self.energy.pressureEnergy + regularization

    def compute_bendingForce(self):
        # self.energy.bendingForce = egrad(
        #     self.compute_meanCurvature)(self.radius_v, self.height_v)
        gradient = grad(self.compute_totalEnergy)
        return gradient
        # print(gradient(np.append(self.radius_v, self.height_v)).shape)
        # self.force.bendingForce = -gradient(self.join(self.radius_v, self.height_v))
        # self.force.bendingForce[0] = 0
        # self.force.bendingForce[1] = 0
        # self.force.bendingForce[10] = 0
        # self.force.bendingForce[11] = 0
        # self.force.bendingForce[12] = 0
        # self.force.bendingForce[13] = 0
        # self.force.bendingForce[22] = 0
        # self.force.bendingForce[23] = 0
        # self.force.bendingForce = egrad(
        #     self.compute_meanCurvature)(self.radius_v)

    def split(self, first_second, partition):
        if len(partition) != 2:
            raise Exception("Can only partition 2!")
        if np.sum(partition) != first_second.shape[0]:
            raise Exception("Partition number does not agree!")
        # twiceN = first_second.shape[0]
        first = first_second[:partition[0]]
        second = first_second[partition[0]:]
        return first, second

    def join(self, first, second):
        return np.append(first, second)
