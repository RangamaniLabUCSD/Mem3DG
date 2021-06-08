// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include <iostream>

#include <gtest/gtest.h>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>

#include "mem3dg/solver/mesh.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class ForceCalculationTest : public testing::Test {
protected:
  // initialize mesh and vpg
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  Eigen::Matrix<size_t, Eigen::Dynamic, 3> topologyMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  Parameters p;
  Options o;
  double h = 0.0001;

  ForceCalculationTest() {

    o.isReducedVolume = true;
    o.isConstantOsmoticPressure = false;
    o.isConstantSurfaceTension = true;
    o.isProteinVariation = true;
    o.isShapeVariation = true;
    o.isFloatVertex = true;
    o.boundaryConditionType = "none";

    o.isEdgeFlip = false;
    o.isSplitEdge = false;
    o.isCollapseEdge = false;
    o.isVertexShift = false;

    p.pt.resize(3, 1);
    p.pt << 0, 0, 1;
    p.protein0.resize(4, 1);
    p.protein0 << 1, 1, 0.7, 0.2;

    p.Kb = 8.22e-5;
    p.Kbc = 0;
    p.H0c = -1;

    p.Ksg = 1e-4;
    p.A_res = 0;
    p.epsilon = -1e-5;

    p.Kv = 0.02;
    p.V_res = 0;
    p.Vt = 0.6;
    p.cam = -1;

    p.Bc = 1;

    p.eta = 0.001;

    p.gamma = 0;
    p.temp = 0;

    p.Kst = 0;
    p.Ksl = 0;
    p.Kse = 0;

    p.Kf = 0;
    p.conc = -1;
    p.height = 0;
    p.radius = -1;
    p.lambdaSG = 0;
    p.lambdaV = 0;

    // Create mesh and geometry objects
    std::tie(topologyMatrix, vertexMatrix) = getIcosphereMatrix(1, 3);
  }
};

/**
 * @brief Test whether passive force is conservative: result need to be the same
 * when computed twice
 *
 */
TEST_F(ForceCalculationTest, ConsistentForcesTest) {
  // Instantiate system object
  size_t nSub = 0;
  mem3dg::System f(topologyMatrix, vertexMatrix, p, o, nSub);
  // First time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX3D vectorForces1 = f.F.toMatrix(f.F.mechanicalForceVec);
  EigenVectorX1D chemicalPotential1 = f.F.toMatrix(f.F.chemicalPotential);
  EigenVectorX3D regularizationForce1 = f.F.toMatrix(f.F.regularizationForce);

  // Second time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX3D vectorForces2 = f.F.toMatrix(f.F.mechanicalForceVec);
  EigenVectorX1D chemicalPotential2 = f.F.toMatrix(f.F.chemicalPotential);
  EigenVectorX3D regularizationForce2 = f.F.toMatrix(f.F.regularizationForce);

  // Comparison of 2 force calculations
  ASSERT_TRUE((vectorForces1 - vectorForces2).norm() < 1e-12);
  ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to decrease in
 * energy
 *
 */
TEST_F(ForceCalculationTest, ConsistentForceEnergy) {

  // initialize the system
  size_t nSub = 0;
  mem3dg::System f(topologyMatrix, vertexMatrix, p, o, nSub);

  // initialize variables
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  const EigenVectorX3D current_pos = f.F.toMatrix(f.vpg->inputVertexPositions);
  const EigenVectorX1D current_proteinDensity = f.F.toMatrix(f.proteinDensity);
  const double tolerance = 1e-4;
  double expectedEnergyDecrease = 0;
  double actualEnergyDecrease = 0;
  double difference = 0;

  f.updateVertexPositions();
  f.computeFreeEnergy();
  Energy previousE{f.E};

  f.computePhysicalForces();

  // bending force
  f.proteinDensity.raw() = current_proteinDensity;
  f.F.toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.F.mask(f.F.toMatrix(f.F.bendingForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * (f.F.mask(f.F.toMatrix(f.F.bendingForceVec)).array() *
           f.F.mask(f.F.toMatrix(f.F.bendingForceVec)).array())
              .sum();
  actualEnergyDecrease = -f.E.BE + previousE.BE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.BE <= previousE.BE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of bending force: " << difference;

  // bending potential
  f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity + h * f.P.Bc * f.F.bendingPotential.raw();
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.P.Bc *
      (f.F.bendingPotential.raw().array() * f.F.bendingPotential.raw().array())
          .sum();
  actualEnergyDecrease = -f.E.BE + previousE.BE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.BE <= previousE.BE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of bending potential: " << difference;

  // capillary force
  f.proteinDensity.raw() = current_proteinDensity;
  f.F.toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.F.mask(f.F.toMatrix(f.F.capillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * (f.F.mask(f.F.toMatrix(f.F.capillaryForceVec)).array() *
           f.F.mask(f.F.toMatrix(f.F.capillaryForceVec)).array())
              .sum();
  actualEnergyDecrease = -f.E.sE + previousE.sE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.sE <= previousE.sE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of capillary force: " << difference;

  // osmotic force
  f.proteinDensity.raw() = current_proteinDensity;
  f.F.toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.F.mask(f.F.toMatrix(f.F.osmoticForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * (f.F.mask(f.F.toMatrix(f.F.osmoticForceVec)).array() *
           f.F.mask(f.F.toMatrix(f.F.osmoticForceVec)).array())
              .sum();
  actualEnergyDecrease = -f.E.pE + previousE.pE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.pE <= previousE.pE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of osmotic force: " << difference;

  // adsorption force
  f.proteinDensity.raw() = current_proteinDensity;
  f.F.toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * (f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec)).array() *
           f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec)).array())
              .sum();
  actualEnergyDecrease = -f.E.aE + previousE.aE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.aE <= previousE.aE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of adsorption force: " << difference;

  // adsorption potential
  f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity + h * f.P.Bc * f.F.adsorptionPotential.raw();
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease = h * (f.F.toMatrix(f.F.adsorptionPotential).array() *
                                f.F.toMatrix(f.F.adsorptionPotential).array())
                                   .sum();
  actualEnergyDecrease = -f.E.aE + previousE.aE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.aE <= previousE.aE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of adsorption potential: " << difference;

  // line tension force
  f.proteinDensity.raw() = current_proteinDensity;
  f.F.toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * (f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec)).array() *
           f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec)).array())
              .sum();
  actualEnergyDecrease = -f.E.dE + previousE.dE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.dE <= previousE.dE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of line tension force: " << difference;

  // diffusion potential
  f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity + h * f.P.Bc * f.F.diffusionPotential.raw();
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease = h * (f.F.toMatrix(f.F.diffusionPotential).array() *
                                f.F.toMatrix(f.F.diffusionPotential).array())
                                   .sum();
  actualEnergyDecrease = -f.E.dE + previousE.dE;
  difference = abs((expectedEnergyDecrease - actualEnergyDecrease) /
                   actualEnergyDecrease);
  ASSERT_TRUE(f.E.dE <= previousE.dE);
  ASSERT_TRUE(difference < tolerance)
      << "Relative error of diffusion potential: " << difference;

}; // namespace ddgsolver
} // namespace mem3dg
