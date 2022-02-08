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
    parameters = Parameters()

    energy = Energy()
    
    force = Force()

    arclength_v = np.array([])

    area_v = np.array([])
    
    volume_v = np.array([])

    turningAngle_v = np.array([])

    radius_v = np.array([])
    radius_e = np.array([])

    height_v = np.array([])

    meanCurvature_v = np.array([])
    radialCurvature_v = np.array([])
    axialCurvature_v = np.array([])

    def __init__(self, radius, height, parameters):
        # self.arclength = arclength
        # self.turning_angle = turning_angle
        self.radius_v = radius
        self.height_v = height
        self.parameters = parameters
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

        self.meanCurvature_v = (self.radialCurvature_v + self.axialCurvature_v) * 0.5

        self.area_v = 2 * np.pi * radius * self.arclength_v
        
        self.volume_v = np.pi * radius_e**2 * dh 

        self.energy.bendingEnergy = self.parameters.Kb * np.sum((self.meanCurvature_v - 10)**2 * self.area_v)
        
        self.energy.surfaceEnergy = self.parameters.Ksg * (np.sum(self.area_v) - self.parameters.At)**2
        
        self.energy.pressureEnergy = self.parameters.Kv * (np.sum(self.volume_v) - self.parameters.Vt)**2

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
