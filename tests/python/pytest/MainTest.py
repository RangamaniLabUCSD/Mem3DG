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

from typing_extensions import assert_type
import pytest
import pymem3dg as dg
import polyscope
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import matplotlib.pyplot as plt
from functools import partial
import pymem3dg.broilerplate as dg_broil
import numpy as np


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

    polyscope.init()
    dg_vis.polyscopeStyle()

    def test_shape_and_protein_variation(self):
        """test simulation with both shape and protein variation"""
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
        p.point.index = 0
        p.external.setForm(
            partial(dg_broil.prescribeGaussianPointForce, Kf=0.005, std=0.02)
        )
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
        """test simulation with only shape variation

        Returns:
            pymem3dg.Parameters: parameter used in integration
        """
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
        p.point.index = 0
        p.external.setForm(
            partial(dg_broil.prescribeGaussianPointForce, Kf=0.005, std=0.02)
        )
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
        """test simulation with only protein variation

        Returns:
            pymem3dg.Parameters: parameter used in integration
        """
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
        """test runs for the convenience function used for visualizing .nc trajectory file"""
        p = self.test_shape_and_protein_variation()
        dg_vis.animate(
            trajNc=self.trajFile,
            parameters=p,
            showBasics=True,
            showForce=True,
            showPotential=True,
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
        """test runs for the visualization of .ply file"""
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

    def test_mesh_generation(self):
        """test mesh generation functions"""
        face, vertex = dg.getTetrahedron()
        face, vertex = dg.getDiamond(dihedral=np.pi / 3)
        face, vertex = dg.getHexagon(radius=1, subdivision=3)
        face, vertex = dg.getCylinder(
            radius=1, radialSubdivision=10, axialSubdivision=10
        )
        face, vertex = dg.getIcosphere(radius=1, subdivision=3)
        face_, vertex_ = dg.linearSubdivide(face=face, vertex=vertex, nSub=2)
        face_, vertex_ = dg.loopSubdivide(face=face, vertex=vertex, nSub=3)

    def test_mesh_reading(self):
        """test mesh reading function"""
        self.test_shape_and_protein_variation()  # generate .ply file for testing
        face, vertex = dg.getFaceAndVertexMatrix(self.plyFile)
        face, vertex = dg.processSoup(self.plyFile)
        elements = dg.getRichDataElementName(self.plyFile)
        properties = dg.getRichDataPropertyName(self.plyFile, "vertex")
        H = dg.getRichData(self.plyFile, "vertex", "meanCurvature")

    def test_mesh_marking(self):
        """test mesh marking functions"""
        face, vertex = dg.getIcosphere(radius=1, subdivision=3)
        faceData = np.zeros(np.shape(face)[0])
        polyscope.remove_all_structures()
        ps_mesh = polyscope.register_surface_mesh("my mesh", vertex, face)
        dg_vis.setPolyscopePermutations(psmesh=ps_mesh, vertex=vertex, face=face)

        def test_locations(embedded_point, accountedCoordinate):
            centerFace, centerBary = dg.getFaceSurfacePointClosestToEmbeddedCoordinate(
                faceMatrix=face,
                vertexMatrix=vertex,
                embeddedCoordinate=embedded_point,
                accountedCoordinate=accountedCoordinate,
            )
            ps_point = polyscope.register_point_cloud(
                f"cloud{embedded_point}", np.array([embedded_point]), radius=0.1
            )
            faceData[centerFace] = 1
            ps_mesh.add_scalar_quantity(
                "center", faceData, defined_on="faces", enabled=True
            )

        test_locations([2, 2, 2], [True, True, True])
        test_locations([2, 2, -2], [True, True, True])
        test_locations([0, 0, -2], [True, True, True])
        test_locations([-1, -4, -2], [True, True, True])
        # polyscope.show()

    def test_parameter_loading(self):
        """test broiler plate example form functions used in parameter loading
        """
        face, vertex = dg.getIcosphere(radius=1, subdivision=3)
        p = dg.Parameters()
        p.point.index = 0
        system = dg.System(face, vertex, p)
        system.initialize()
        system.computeGeodesicDistance()
        system.parameters.protein.setForm(
            partial(
                dg_broil.prescribeGeodesicPoteinDensityDistribution,
                sharpness=20,
                radius=0.1,
            )
        )
        system.prescribeProteinDensityDistribution()
        system.parameters.external.setForm(
            partial(dg_broil.prescribeGaussianPointForce, Kf=0.01, std=1)
        )
        system.prescribeExternalForce()
        system.parameters.external.setForm(
            partial(dg_broil.prescribePeriodicForceOnCylinder, Kf=0.01, freq=10)
        )
        system.prescribeExternalForce()
