import autograd.numpy as np
from autograd import elementwise_grad as egrad
from autograd import grad
from ._core import *


class Parameters:
    Kb = 0
    Kbc = 0
    Ksg = 0
    At = 0
    epsilon = 0
    Kv = 0

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


class Force:
    bendingForce = np.array([])
    surfaceForce = np.array([])
    pressureForce = np.array([])
    adsorptionForce = np.array([])

    def __init__(self) -> None:
        pass


class System:
    parameter = Parameters()

    energy = Energy()

    arclength_v = np.array([])

    turningAngle_v = np.array([])

    radius_v = np.array([])
    radius_e = np.array([])

    height_v = np.array([])

    meanCurvature_v = np.array([])
    radialCurvature_v = np.array([])
    axialCurvature_v = np.array([])

    def __init__(self, radius, height, parameter):
        # self.arclength = arclength
        # self.turning_angle = turning_angle
        self.radius_v = radius
        self.height_v = height
        self.parameter = parameter
        self.compute_totalEnergy(self.join(self.radius_v, self.height_v))

    def split(self, radius_height):
        twiceN = radius_height.shape[0]
        radius = radius_height[:int(twiceN/2)]
        height = radius_height[int(twiceN/2):]
        return radius, height

    def join(self, radius, height):
        return np.append(radius, height)

    def compute_totalEnergy(self, radius_height):
        radius, height = self.split(radius_height)
        n = radius.shape[0]

        dr = np.diff(radius)
        dh = np.diff(height)
        ds = np.sqrt(dr**2 + dh**2)
        s1 = np.append(ds, 0.0)
        s2 = np.append(0.0, ds)
        self.arclength_v = (s1 + s2) * 0.5

        psi = np.arctan2(dh, dr)

        phi = np.diff(psi)
        phi = np.append(0.0, phi)
        phi = np.append(phi, 0.0)
        self.turningAngle_v = phi

        self.axialCurvature_v = self.turningAngle_v / self.arclength_v

        radius_e = (radius[0:n-1] + radius[1:n]) * 0.5
        radialCurvature_e = np.sin(psi) / radius_e * ds
        rc1 = np.append(radialCurvature_e, 0.0)
        rc2 = np.append(0.0, radialCurvature_e)
        self.radialCurvature_v = (rc1 + rc2) * 0.5 / self.arclength_v

        self.meanCurvature_v = self.radialCurvature_v + self.axialCurvature_v

        self.energy.bendingEnergy = np.sum(self.meanCurvature_v * self.meanCurvature_v * self.arclength_v)

        return self.energy.bendingEnergy

    def compute_bendingForce(self):
        # self.energy.bendingForce = egrad(
        #     self.compute_meanCurvature)(self.radius_v, self.height_v)
        gradient = grad(self.compute_totalEnergy)
        # print(gradient(np.append(self.radius_v, self.height_v)).shape)
        self.energy.bendingForce = gradient(self.join(self.radius_v, self.height_v))
        # self.energy.bendingForce = egrad(
        #     self.compute_meanCurvature)(self.radius_v)
