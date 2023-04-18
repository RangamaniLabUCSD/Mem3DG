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
import pymem3dg as dg
import pymem3dg.read.netcdf as dg_nc
import pymem3dg.read.meshop as dg_meshop
import pymem3dg.read.mesh as dg_mesh
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import pymem3dg.broilerplate as dg_broil

from cmath import sqrt
import polyscope
import matplotlib.pyplot as plt
from functools import partial
import numpy as np
from pathlib import Path


class TestLoadCheck(object):
    def test_docs(self):
        print(dg.__doc__)
        assert True


class TestInitialization(object):
    face, vertex = dg_meshop.getIcosphere(1, 3)
    vertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1)
    geometry = dg.Geometry(face, vertex, vertex)
    proteinDensity = np.ones(np.shape(vertex)[0]) * 0.1
    velocity = np.zeros(np.shape(vertex))
    initialConditions = {
        "geometry": geometry,
        "proteinDensity": proteinDensity,
        "velocity": velocity,
    }
    notableVertex = [False] * np.shape(vertex)[0]
    notableVertex[10] = True
    notableVertex[50] = True

    def test_geometry(self):
        # test matrix based initialization
        g1 = dg.Geometry(self.face, self.vertex, self.vertex, self.notableVertex)
        g2 = dg.Geometry(self.face, self.vertex, self.vertex)
        # test the notable vertex argument
        assert (g1.getNotableVertex() == self.notableVertex).all()
        assert g2.getNotableVertex()[0]

    def test_system(self):
        p = dg.Parameters()
        p.variation.isShapeVariation = True
        p.bending.Kbc = 0.1
        p.bending.H0c = 10
        p.tension.form = partial(dg_broil.constantSurfaceTensionModel, tension=0.5)
        p.osmotic.form = partial(dg_broil.constantOsmoticPressureModel, pressure=0.01)
        arguments = self.initialConditions
        arguments["parameters"] = p
        s1 = dg.System(**arguments)
        assert (
            s1.getGeometry().getVertexMatrix() == self.geometry.getVertexMatrix()
        ).all()

    def test_parameter_loading(self):
        """test broiler plate example form functions used in parameter loading"""
        face, vertex = dg_meshop.getIcosphere(radius=1, subdivision=3)
        p = dg.Parameters()
        geometry = dg.Geometry(face, vertex, vertex, self.notableVertex)
        system = dg.System(geometry, p)
        system.initialize()

        # protein
        system.parameters.protein.prescribeProteinDensityDistribution = partial(
            dg_broil.prescribeGeodesicPoteinDensityDistribution,
            sharpness=20,
            radius=0.5,
        )
        system.initialize()

        ## uncomment the snippet to visualize
        # polyscope.init()
        # psmesh = dg_vis.visualizeGeometry(geometry)
        # psmesh.add_scalar_quantity("proteindensity", system.getProteinDensity(), cmap = "coolwarm")
        # polyscope.show()

        # external force
        system.parameters.external.form = partial(
            dg_broil.prescribeGaussianPointForce, Kf=0.01, std=1, tau=100
        )
        system.prescribeExternalForce()
        system.parameters.external.form = partial(
            dg_broil.prescribePeriodicForceOnCylinder, Kf=0.01, freq=10
        )
        system.prescribeExternalForce()

        # osmotic pressure
        system.parameters.osmotic.form = partial(
            dg_broil.constantOsmoticPressureModel, pressure=0.01
        )
        system.initialize()
        assert system.getForces().getOsmoticPressure() == 0.01
        assert (
            system.getEnergy().pressureEnergy
            == -0.01 * system.getGeometry().getVolume()
        )

        # surface tension
        system.parameters.tension.form = partial(
            dg_broil.constantSurfaceTensionModel, tension=0.01
        )
        system.initialize()
        assert system.getForces().getSurfaceTension() == 0.01
        assert (
            system.getEnergy().surfaceEnergy
            == 0.01 * system.getGeometry().getSurfaceArea()
        )

        # test read
        assert system.parameters.osmotic.form(2) == (0.01, -0.01 * 2)


class TestExampleIntegration(object):
    outputDir = str(Path(__file__).parent.absolute())
    trajFile = outputDir + "/traj.nc"
    plyFile = outputDir + "/test.ply"
    initialization = TestInitialization()

    def test_shape_and_protein_variation(self):
        """test simulation with both shape and protein variation"""
        p = dg.Parameters()
        p.variation.isProteinVariation = True
        p.variation.isShapeVariation = True
        p.bending.Kbc = 0.1
        p.bending.H0c = 10
        p.tension.form = partial(dg_broil.constantSurfaceTensionModel, tension=0.5)
        p.osmotic.form = partial(dg_broil.constantOsmoticPressureModel, pressure=0.01)
        p.dirichlet.eta = p.bending.Kb
        p.proteinMobility = 1
        p.spring.Kst = 1
        p.external.form = partial(
            dg_broil.prescribeGaussianPointForce, Kf=0.005, std=0.02, tau=100
        )
        arguments = self.initialization.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(ifMutateMesh=False, ifMute=False)
        g.computeTotalEnergy()
        g.computeConservativeForcing()
        g.addNonconservativeForcing()
        g.saveRichData(self.plyFile)

        fe = dg.Euler(
            system=g,
            characteristicTimeStep=1e-3,
            totalTime=1e-1,
            savePeriod=1e-2,
            tolerance=1e-10,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        success = fe.integrate()
        assert success == True
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
        p.tension.form = partial(dg_broil.constantSurfaceTensionModel, tension=0.5)
        p.osmotic.form = partial(dg_broil.constantOsmoticPressureModel, pressure=0.01)
        p.spring.Kst = 1
        p.external.form = partial(
            dg_broil.prescribeGaussianPointForce, Kf=0.005, std=0.02, tau=100
        )
        arguments = self.initialization.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(ifMutateMesh=False, ifMute=False)
        fe = dg.Euler(
            system=g,
            characteristicTimeStep=1e-3,
            totalTime=1e-1,
            savePeriod=1e-2,
            tolerance=1e-10,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        success = fe.integrate()
        assert success == True
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
        arguments = self.initialization.initialConditions
        arguments["parameters"] = p
        g = dg.System(**arguments)
        g.initialize(ifMutateMesh=False, ifMute=False)
        fe = dg.Euler(
            system=g,
            characteristicTimeStep=0.1,
            totalTime=1,
            savePeriod=0.1,
            tolerance=1e-10,
            outputDirectory=self.outputDir,
        )
        fe.ifPrintToConsole = True
        fe.ifOutputTrajFile = True
        success = fe.integrate()
        assert success == True
        return p


class TestContinuation(object):
    trajectory = TestExampleIntegration()
    parameters = trajectory.test_shape_and_protein_variation()

    def test_initialization(self):
        frame = 2
        face, vertex = dg_nc.getFaceAndVertexMatrix(
            trajNc=self.trajectory.trajFile, frame=frame
        )
        g2 = dg.Geometry(face, vertex, vertex)
        g3 = dg.Geometry(self.trajectory.trajFile, frame)
        assert g3.getSurfaceArea() == g2.getSurfaceArea()
        g = dg.Geometry(self.trajectory.trajFile, 0)
        s2 = dg.System(g, self.trajectory.trajFile, 0, self.parameters)


class TestMeshIO(object):
    def test_mesh_generation(self):
        """test mesh generation functions"""
        face, vertex = dg_meshop.getTetrahedron()
        face, vertex = dg_meshop.getDiamond(dihedral=np.pi / 3)
        face, vertex = dg_meshop.getHexagon(radius=1, subdivision=3)
        face, vertex = dg_meshop.getCylinder(
            radius=1, radialSubdivision=10, axialSubdivision=10
        )
        face, vertex = dg_meshop.getIcosphere(radius=1, subdivision=3)
        face_, vertex_ = dg_meshop.linearSubdivide(face=face, vertex=vertex, nSub=2)
        face_, vertex_ = dg_meshop.loopSubdivide(face=face, vertex=vertex, nSub=3)

    def test_mesh_reading(self):
        """test mesh reading function"""
        g = TestExampleIntegration()
        g.test_shape_and_protein_variation()  # files for testing

        face, vertex = dg_mesh.getFaceAndVertexMatrix(g.plyFile)
        face, vertex = dg_mesh.processSoup(g.plyFile)
        elements = dg_mesh.getDataElementName(g.plyFile)
        properties = dg_mesh.getDataPropertyName(g.plyFile, "vertex")
        H = dg_mesh.getData(g.plyFile, "vertex", "meanCurvature")

        face, vertex = dg_nc.getFaceAndVertexMatrix(trajNc=g.trajFile, frame=3)
        vel = dg_nc.getData(
            trajNc=g.trajFile,
            frame=3,
            group="Trajectory",
            variable="velocities",
            num_col=3,
        )
        phi = dg_nc.getData(
            trajNc=g.trajFile,
            frame=3,
            group="Trajectory",
            variable="proteindensity",
            num_col=1,
        )

        # dg_nc.crop(trajnc=g.trajFile, limit=(0, 20))

    def test_mesh_marking_1(self):
        """test mesh marking functions"""
        face, vertex = dg_meshop.getIcosphere(radius=1, subdivision=3)
        faceData = np.zeros(np.shape(face)[0])

        # polyscope.init()
        # polyscope.remove_all_structures()
        # ps_mesh = polyscope.register_surface_mesh("my mesh", vertex, face)
        # dg_vis.setPolyscopePermutations(psmesh=ps_mesh, vertex=vertex, face=face)
        def test_locations(embedded_point, accountedCoordinate):
            (
                centerFace,
                centerBary,
            ) = dg_meshop.getFaceSurfacePointClosestToEmbeddedCoordinate(
                faceMatrix=face,
                vertexMatrix=vertex,
                embeddedCoordinate=embedded_point,
                accountedCoordinate=accountedCoordinate,
            )
            center = dg_meshop.getVertexClosestToEmbeddedCoordinate(
                vertexMatrix=vertex,
                embeddedCoordinate=embedded_point,
                accountedCoordinate=accountedCoordinate,
            )
            # ps_point = polyscope.register_point_cloud(
            #     f"cloud{embedded_point}", np.array([embedded_point]), radius=0.1
            # )
            # faceData[centerFace] = 1
            # ps_mesh.add_scalar_quantity(
            #     "center", faceData, defined_on="faces", enabled=True
            # )
            return centerFace, center

        test_locations([2, 2, 2], [True, True, True])
        test_locations([2, 2, -2], [True, True, True])
        test_locations([0, 0, -2], [True, True, True])
        test_locations([-1, -4, -2], [True, True, True])

    def test_mesh_marking_2(self):
        # test center prescription on open boundary mesh
        face, vertex = dg_meshop.getHexagon(radius=1, subdivision=4)
        vertex = vertex + np.ones(np.shape(vertex))
        notableVertex = dg_meshop.getVertexFurthestFromBoundary(face, vertex)
        assert np.linalg.norm(vertex[notableVertex]) == sqrt(3)


class TestVisualization(object):
    # polyscope.init()
    dg_vis.polyscopeStyle()
    initialization = TestInitialization()
    trajectory = TestExampleIntegration()
    parameters = trajectory.test_shape_and_protein_variation()

    def test_nc_visual(self):
        """test runs for the convenience function used for visualizing .nc trajectory file"""
        # with parameters
        dg_vis.animate(
            trajNc=self.trajectory.trajFile,
            parameters=self.parameters,
            showBasics=True,
            showForce=True,
            showPotential=True,
        )
        # without parameters
        dg_vis.animate(
            trajNc=self.trajectory.trajFile,
            showBasics=True,
            showForce=True,
            showPotential=True,
        )
        _, ax = plt.subplots(4)
        dg_vis.plotProteinDensity(ax[0], self.trajectory.trajFile)
        dg_vis.plotEnergy(
            ax=ax[1],
            trajFile=self.trajectory.trajFile,
            parameters=self.parameters,
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
            self.trajectory.trajFile,
            self.parameters,
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
            self.trajectory.trajFile,
            self.parameters,
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
        dg_vis.visualizePly(
            self.trajectory.plyFile,
            "proteinDensity",
            "forceMask",
            "proteinMask",
            "meanCurvature",
            "notableVertex",
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

    def test_geometry_visual(self):
        dg_vis.visualizeGeometry(self.initialization.geometry)
        # polyscope.show()
