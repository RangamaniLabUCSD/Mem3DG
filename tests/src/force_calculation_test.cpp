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

#include "mem3dg/mesh_io.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace solver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class ForceCalculationTest : public testing::Test {
protected:
  // initialize mesh and vpg
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> topologyMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  Parameters p;
  Options o;
  double h = 0.1;

  ForceCalculationTest() {
    o.isFloatVertex = true;

    // o.isEdgeFlip = false;
    // o.isSplitEdge = false;
    // o.isCollapseEdge = false;
    o.isVertexShift = false;

    p.variation.isShapeVariation = true;
    p.variation.isProteinVariation = true;
    p.pt.resize(3, 1);
    p.pt << 0, 0, 1;
    p.protein0.resize(4, 1);
    p.protein0 << 1, 1, 0.7, 0.2;
    p.sharpness = 3;

    p.bending.Kb = 8.22e-5;
    p.bending.Kbc = 0;
    p.bending.H0c = -1;

    p.tension.isConstantSurfaceTension = true;
    p.tension.Ksg = 1e-2;
    p.tension.A_res = 0;

    p.adsorption.epsilon = -1e-2;

    p.osmotic.isPreferredVolume = false;
    p.osmotic.isConstantOsmoticPressure = true;
    p.osmotic.Kv = 1e-2;
    p.osmotic.V_res = 0;
    p.osmotic.Vt = -1;
    p.osmotic.cam = -1;
    p.osmotic.n = 1;

    p.boundary.shapeBoundaryCondition = "roller";
    p.boundary.proteinBoundaryCondition = "pin";

    p.Bc = 1;

    p.dirichlet.eta = 0.001;

    p.dpd.gamma = 0;

    p.temp = 0;

    p.Kst = 0;
    p.Ksl = 0;
    p.Kse = 0;

    p.external.Kf = 0;
    p.external.conc = -1;
    p.external.height = 0;
    p.radius = -1;
    p.lambdaSG = 0;
    p.lambdaV = 0;

    // Create mesh and geometry objects
    std::tie(topologyMatrix, vertexMatrix) = getCylinderMatrix(1, 10, 10);
  }
};

/**
 * @brief Test whether passive force is conservative: result need to be the same
 * when computed twice
 *
 */
TEST_F(ForceCalculationTest, ConsistentForcesTest) {
  // Instantiate system object
  std::size_t nSub = 0;
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, p, o, nSub);
  // First time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX3dr mechanicalForceVec1 = toMatrix(f.forces.mechanicalForceVec);
  EigenVectorX1d chemicalPotential1 = toMatrix(f.forces.chemicalPotential);
  EigenVectorX3dr regularizationForce1 = toMatrix(f.forces.regularizationForce);

  // Second time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX3dr mechanicalForceVec2 = toMatrix(f.forces.mechanicalForceVec);
  EigenVectorX1d chemicalPotential2 = toMatrix(f.forces.chemicalPotential);
  EigenVectorX3dr regularizationForce2 = toMatrix(f.forces.regularizationForce);

  // Comparison of 2 force calculations
  EXPECT_TRUE(mechanicalForceVec1.isApprox(mechanicalForceVec2));
  EXPECT_TRUE(chemicalPotential1.isApprox(chemicalPotential2));
  EXPECT_TRUE(regularizationForce1.isApprox(regularizationForce2));
  //   EXPECT_TRUE((mechanicalForceVec1 - mechanicalForceVec2).norm() < 1e-12);
  //   EXPECT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  //   EXPECT_TRUE((regularizationForce1 - regularizationForce2).norm() <
  //   1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to decrease in
 * energy
 *
 */
TEST_F(ForceCalculationTest, ConsistentForceEnergy) {

  // initialize the system
  std::size_t nSub = 0;
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, p, o, nSub);

  // initialize variables
  auto vel_e = gc::EigenMap<double, 3>(f.velocity);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  const EigenVectorX3dr current_pos = toMatrix(f.vpg->inputVertexPositions);
  const EigenVectorX1d current_proteinDensity = toMatrix(f.proteinDensity);
  const double tolerance = 0.05;
  double expectedEnergyDecrease = 0;
  double actualEnergyDecrease = 0;
  double difference_h = 0;
  double difference_xh = 0;
  std::size_t stepFold = 2;
  double expectRate = 2;

  f.updateVertexPositions();
  f.computeFreeEnergy();
  Energy previousE{f.energy};

  f.computePhysicalForces();

  // bending force
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.BE + previousE.BE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.BE <= previousE.BE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "bending force: (expected - actual) / expected = "
      << difference_h / abs(actualEnergyDecrease);

  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + stepFold * h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.BE + previousE.BE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "bending force: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // bending potential
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity + h * f.forces.maskProtein(f.forces.bendingPotential.raw());
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.bendingPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.BE + previousE.BE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.BE <= previousE.BE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "bending potential: (expected - actual) / expected: "
      << difference_h / abs(actualEnergyDecrease);

  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h * f.forces.maskProtein(f.forces.bendingPotential.raw());
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h * f.forces.maskProtein(f.forces.bendingPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.BE + previousE.BE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "bending potential: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // capillary force
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.sE + previousE.sE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.sE <= previousE.sE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "capillary force: (expected - actual) / expected = "
      << difference_h / abs(actualEnergyDecrease);

  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.capillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.sE + previousE.sE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "capillary force: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // osmotic force
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.pE + previousE.pE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.pE <= previousE.pE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "osmotic force: (expected - actual) / expected = "
      << difference_h / abs(actualEnergyDecrease);

  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + stepFold * h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.pE + previousE.pE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "osmotic force: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // adsorption force
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.aE + previousE.aE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.aE <= previousE.aE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "adsorption force: (expected - actual) / expected = "
      << difference_h / abs(actualEnergyDecrease);

  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.aE + previousE.aE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "adsorption force: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // adsorption potential
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() = current_proteinDensity +
                           h * f.forces.maskProtein(f.forces.adsorptionPotential.raw());
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.adsorptionPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.aE + previousE.aE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.aE <= previousE.aE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "adsorption potential: (expected - actual) / expected = "
      << difference_h / abs(actualEnergyDecrease);

  //   toMatrix(f.vpg->inputVertexPositions) = current_pos;
  //   f.proteinDensity.raw() =
  //       current_proteinDensity +
  //       stepFold * h * f.F.maskProtein(f.F.adsorptionPotential.raw());
  //   f.updateVertexPositions(false);
  //   f.computeFreeEnergy();
  //   expectedEnergyDecrease =
  //       stepFold * h *
  //       f.F.maskProtein(f.F.adsorptionPotential.raw()).squaredNorm();
  //   actualEnergyDecrease = -f.E.aE + previousE.aE;
  //   difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  //   EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
  //               tolerance)
  //       << "adsorption potential: difference_xh / difference_h = "
  //       << difference_xh / difference_h;

  // line tension force
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.dE + previousE.dE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.dE <= previousE.dE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "line tension force = (expected - actual) / expected: "
      << difference_h / abs(actualEnergyDecrease);

  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec));
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.dE + previousE.dE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "line tension force: difference_xh / difference_h = "
      << difference_xh / difference_h;

  // diffusion potential
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() = current_proteinDensity +
                           h * f.forces.maskProtein(f.forces.diffusionPotential.raw());
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.diffusionPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.dE + previousE.dE;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.dE <= previousE.dE);
  EXPECT_TRUE(difference_h < tolerance * abs(actualEnergyDecrease))
      << "diffusion potential = (expected - actual) / expected: "
      << difference_h / abs(actualEnergyDecrease);

  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h * f.forces.maskProtein(f.forces.diffusionPotential.raw());
  f.updateVertexPositions(false);
  f.computeFreeEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskProtein(f.forces.diffusionPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.dE + previousE.dE;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "diffusion potential: difference_xh / difference_h = "
      << difference_xh / difference_h;
};
} // namespace solver
} // namespace mem3dg
