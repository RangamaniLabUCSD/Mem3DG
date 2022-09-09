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

from numpy import polyfit
import pytest
import pymem3dg as dg
import polyscope
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import matplotlib.pyplot as plt


class TestLoadCheck(object):
    def test_docs(self):
        print(dg.__doc__)
        assert 1 == 1


class TestExampleIntegration(object):
    from pathlib import Path

    outputDir = str(Path(__file__).parent.absolute())
    trajFile = outputDir + "/traj.nc"
    plyFile = outputDir + "/test.ply"
    face, vertex = dg.getIcosphere(1, 3)

    vertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1)
    import numpy as np

    proteinDensity = np.ones(np.shape(vertex)[0]) * 0.1
    velocity = np.zeros(np.shape(vertex))
    initialConditions = {
        "topologyMatrix": face,
        "vertexMatrix": vertex,
        "referenceVertexMatrix": vertex,
        "proteinDensity": proteinDensity,
        "velocity": velocity,
    }

    def test_shape_and_protein_variation(self):
        def point(vertexPositions, vertexDualAreas, time, geodesicDistances):
            direction = dg_util.rowwiseNormalize(vertexPositions)
            magnitude = 0.005 * dg_util.gaussianDistribution(geodesicDistances, 0, 0.02)
            return dg_util.rowwiseScaling(magnitude, direction)

        p = dg.Parameters()
        p.variation.isProteinVariation = True
        p.variation.isShapeVariation = True
        p.bending.Kbc = 0.1
        p.bending.H0c = 10
        p.tension.isConstantSurfaceTension = True
        p.tension.Ksg = 0.5
        p.osmotic.isConstantOsmoticPressure = True
        p.osmotic.Kv = 0.01
        p.dirichlet.eta = p.bending.Kb
        p.proteinMobility = 1
        p.spring.Kst = 1
        p.point.pt = [0, 0, 10]
        p.external.setForm(point)
        arguments = self.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(nMutation=0, ifMute=False)
        g.computeTotalEnergy()
        g.computeConservativeForcing()
        g.addNonconservativeForcing()
        g.saveRichData(self.plyFile)

        fe = dg.Euler(
            system=g,
            characteristicTimeStep=1e-3,
            totalTime=1,
            savePeriod=1e-1,
            tolerance=1e-6,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        sucess = fe.integrate()
        assert sucess == True
        return p

    def test_shape_variation(self):
        def point(vertexPositions, vertexDualAreas, time, geodesicDistances):
            direction = dg_util.rowwiseNormalize(vertexPositions)
            magnitude = 0.005 * dg_util.gaussianDistribution(geodesicDistances, 0, 0.02)
            return dg_util.rowwiseScaling(magnitude, direction)

        p = dg.Parameters()
        p.variation.isShapeVariation = True
        p.variation.isProteinVariation = False
        p.bending.Kbc = 0.1
        p.bending.H0c = 10
        p.tension.isConstantSurfaceTension = True
        p.tension.Ksg = 0.5
        p.osmotic.isConstantOsmoticPressure = True
        p.osmotic.Kv = 0.01
        p.spring.Kst = 1
        p.point.pt = [0, 0, 10]
        p.external.setForm(point)
        arguments = self.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(nMutation=0, ifMute=False)
        fe = dg.Euler(
            system=g,
            characteristicTimeStep=1e-3,
            totalTime=1,
            savePeriod=1e-1,
            tolerance=1e-6,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        sucess = fe.integrate()
        assert sucess == True
        return p

    def test_protein_variation(self):
        p = dg.Parameters()
        p.variation.isShapeVariation = False
        p.variation.isProteinVariation = True
        p.bending.Kb = 8.22e-5
        p.bending.H0c = 10
        p.dirichlet.eta = p.bending.Kb
        p.adsorption.epsilon = -1e-3
        p.proteinMobility = 1
        arguments = self.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(nMutation=0, ifMute=False)
        fe = dg.Euler(
            system=g,
            characteristicTimeStep=0.1,
            totalTime=1,
            savePeriod=0.1,
            tolerance=1e-6,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        sucess = fe.integrate()
        assert sucess == True
        return p

    def test_nc_visual(self):
        p = self.test_shape_and_protein_variation()
        dg_vis.animate(
            trajNc=self.trajFile,
            parameters=p,
            frames=[0, 4, 5]
        )
        _, ax = plt.subplots(4)
        dg_vis.plotProteinDensity(ax[0], self.trajFile, p)
        dg_vis.plotEnergy(
            ax=ax[1],
            trajFile=self.trajFile,
            parameters=p,
            zeroing=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            spontaneousCurvatureEnergy=True,
            externalWork=True,
            deviatoricCurvatureEnergy=True,
            surfaceEnergy=True,
            pressureEnergy=True,
            adsorptionEnergy=True,
            aggregationEnergy=True,
            entropyEnergy=True,
            edgeSpringEnergy=True,
            faceSpringEnergy=True,
            lcrSpringEnergy=True,
            dirichletEnergy=True,
        )
        dg_vis.plotChemicalPotentials(
            ax[2],
            self.trajFile,
            p,
            logScale=True,
            spontaneousCurvaturePotential=True,
            deviatoricCurvaturePotential=True,
            aggregationPotential=True,
            entropyPotential=True,
            dirichletPotential=True,
            adsorptionPotential=True,
        )
        dg_vis.plotMechanicalForces(
            ax[3],
            self.trajFile,
            p,
            logScale=True,
            spontaneousCurvatureForce=True,
            capillaryForce=True,
            osmoticForce=True,
            adsorptionForce=True,
            lineCapillaryForce=True,
            aggregationForce=True,
            externalForce=True,
            entropyForce=True,
            springForce=True,
        )
        # polyscope.show()
        # plt.show()

    def test_ply_visual(self):
        self.test_shape_and_protein_variation()
        dg_vis.visualizePly(
            self.plyFile,
            "proteinDensity",
            "forceMask",
            "proteinMask",
            "meanCurvature",
            "center",
            "gaussianCurvature",
            "spontaneousCurvature",
            "dualArea",
            "spontaneousCurvatureForce",
            "deviatoricCurvatureForce",
            "areaDifferenceForce",
            "capillaryForce",
            "lineCapillaryForce",
            "osmoticForce",
            "adsorptionForce",
            "aggregationForce",
            "externalForce",
            "selfAvoidanceForce",
            "mechanicalForce",
            "dirichletPotential",
            "spontaneousCurvaturePotential",
            "deviatoricCurvaturePotential",
            "adsorptionPotential",
            "aggregationPotential",
            "chemicalPotential",
        )
        # polyscope.show()
