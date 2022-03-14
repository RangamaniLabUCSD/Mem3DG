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

class ForceTest : public testing::Test {
protected:
  // initialize mesh and vpg
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  // std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> topologyMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  Parameters p;
  double h = 0.1;

  ForceTest() {

    p.variation.isShapeVariation = true;
    p.variation.isProteinVariation = true;
    p.variation.radius = -1;
    p.point.isFloatVertex = false;
    p.point.pt.resize(3, 1);
    p.point.pt << 0, 0, 1;
    p.proteinDistribution.protein0.resize(4, 1);
    p.proteinDistribution.profile = "tanh";
    p.proteinDistribution.protein0 << 1, 1, 0.7, 0.2;
    p.proteinDistribution.tanhSharpness = 3;

    p.bending.Kd = 8.22e-5;
    p.bending.Kdc = 8.22e-5;
    p.bending.Kb = 8.22e-5;
    p.bending.Kbc = 0;
    p.bending.H0c = -1;

    p.tension.isConstantSurfaceTension = true;
    p.tension.Ksg = 1e-2;
    p.tension.A_res = 0;
    p.tension.lambdaSG = 0;

    p.adsorption.epsilon = -1e-2;

    p.aggregation.chi = -1e-2;

    p.osmotic.isPreferredVolume = false;
    p.osmotic.isConstantOsmoticPressure = true;
    p.osmotic.Kv = 1e-2;
    p.osmotic.V_res = 0;
    p.osmotic.Vt = -1;
    p.osmotic.cam = -1;
    p.osmotic.n = 1;
    p.osmotic.lambdaV = 0;

    p.boundary.shapeBoundaryCondition = "pin";
    p.boundary.proteinBoundaryCondition = "pin";

    p.proteinMobility = 1;

    p.dirichlet.eta = 0.001;

    p.selfAvoidance.mu = 1e-5;

    p.dpd.gamma = 0;

    p.temperature = 0;

    p.external.Kf = 0;

    // Create mesh and geometry objects
    std::tie(topologyMatrix, vertexMatrix) =
        getCylinderMatrix(1, 10, 10, 5, 0.3);
  }
};

/**
 * @brief Test whether force is conservative: result need to be the same
 * when computed twice
 *
 */
TEST_F(ForceTest, ConservativeForcesTest) {
  // Instantiate system object
  std::size_t nSub = 0;
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, p, nSub);
  // First time calculation of force
  f.computePhysicalForcing();
  f.computeRegularizationForce();
  EigenVectorX3dr mechanicalForceVec1 = toMatrix(f.forces.mechanicalForceVec);
  EigenVectorX1d chemicalPotential1 = toMatrix(f.forces.chemicalPotential);
  // Second time calculation of force
  f.computePhysicalForcing();
  f.computeRegularizationForce();
  EigenVectorX3dr mechanicalForceVec2 = toMatrix(f.forces.mechanicalForceVec);
  EigenVectorX1d chemicalPotential2 = toMatrix(f.forces.chemicalPotential);

  // Comparison of 2 force calculations
  EXPECT_TRUE(mechanicalForceVec1.isApprox(mechanicalForceVec2));
  EXPECT_TRUE(chemicalPotential1.isApprox(chemicalPotential2));
  //   EXPECT_TRUE((mechanicalForceVec1 - mechanicalForceVec2).norm() < 1e-12);
  //   EXPECT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  //   EXPECT_TRUE((regularizationForce1 - regularizationForce2).norm() <
  //   1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to
 * 1. decrease in energy
 * 2. decrease in second order (or exact)
 */
TEST_F(ForceTest, ConsistentForceEnergy) {

  // initialize the system
  std::size_t nSub = 0;
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, p, nSub);

  // initialize variables
  auto vel_e = gc::EigenMap<double, 3>(f.velocity);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  const EigenVectorX3dr current_pos = toMatrix(f.vpg->inputVertexPositions);
  const EigenVectorX1d current_proteinDensity = toMatrix(f.proteinDensity);
  const double tolerance = 0.03;
  double expectedEnergyDecrease = 0;
  double actualEnergyDecrease = 0;
  double difference_h = 0;
  double difference_xh = 0;
  std::size_t stepFold = 2;
  double expectRate = 2;

  f.updateConfigurations();
  f.computeTotalEnergy();
  Energy previousE{f.energy};

  f.computePhysicalForcing();

#pragma region force
  // ==========================================================
  // ================        Bending         ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec));
  f.updateConfigurations(false);
  f.computeBendingEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.bendingEnergy + previousE.bendingEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.bendingEnergy <= previousE.bendingEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "bending force: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.bendingForceVec));
  f.updateConfigurations(false);
  f.computeBendingEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.bendingForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.bendingEnergy + previousE.bendingEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "bending force convergence rate" << std::endl;

  // ==========================================================
  // ================      Deviatoric        ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.deviatoricForceVec));
  f.updateConfigurations(false);
  f.computeDeviatoricEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskForce(toMatrix(f.forces.deviatoricForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.deviatoricEnergy + previousE.deviatoricEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.deviatoricEnergy <= previousE.deviatoricEnergy)
      << "deviatoricPre vs now: " << previousE.deviatoricEnergy
      << f.energy.deviatoricEnergy << "\n";
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "deviatoric force: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;
  //   << f.forces.maskForce(toMatrix(f.forces.deviatoricForceVec)) <<
  //   std::endl;

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.deviatoricForceVec));
  f.updateConfigurations(false);
  f.computeDeviatoricEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.deviatoricForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.deviatoricEnergy + previousE.deviatoricEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(difference_xh < tolerance * (abs(actualEnergyDecrease)))
      << "deviatoric force 1/2: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "deviatoric force convergence rate" << std::endl;

  // ==========================================================
  // ================      Capillary         ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec));
  f.updateConfigurations(false);
  f.computeSurfaceEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskForce(toMatrix(f.forces.capillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.surfaceEnergy + previousE.surfaceEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.surfaceEnergy <= previousE.surfaceEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "capillary force: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec));
  f.updateConfigurations(false);
  f.computeSurfaceEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.capillaryForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.surfaceEnergy + previousE.surfaceEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "capillary force convergence rate" << std::endl;

  // ==========================================================
  // ================      Osmotic           ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos + h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec));
  f.updateConfigurations(false);
  f.computePressureEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.pressureEnergy + previousE.pressureEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.pressureEnergy <= previousE.pressureEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "osmotic force: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec));
  f.updateConfigurations(false);
  f.computePressureEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.osmoticForceVec)).squaredNorm();
  actualEnergyDecrease = -f.energy.pressureEnergy + previousE.pressureEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "osmotic force convergence rate" << std::endl;

  // ==========================================================
  // ================      Adsorption        ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec));
  f.updateConfigurations(false);
  f.computeAdsorptionEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.adsorptionEnergy + previousE.adsorptionEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.adsorptionEnergy <= previousE.adsorptionEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease) + 1e-12))
      << "adsorption force: (expected - actual) / expected = "
      << difference_h / (abs(actualEnergyDecrease));

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec));
  f.updateConfigurations(false);
  f.computeAdsorptionEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.adsorptionEnergy + previousE.adsorptionEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "adsorption force convergence rate" << std::endl;

  // ==========================================================
  // ================      Aggregation       ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.aggregationForceVec));
  f.updateConfigurations(false);
  f.computeAggregationEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskForce(toMatrix(f.forces.aggregationForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.aggregationEnergy + previousE.aggregationEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.aggregationEnergy <= previousE.aggregationEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease) + 1e-12))
      << "aggregation force: (expected - actual) / expected = "
      << difference_h / (abs(actualEnergyDecrease));

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h * f.forces.maskForce(toMatrix(f.forces.aggregationForceVec));
  f.updateConfigurations(false);
  f.computeAggregationEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.aggregationForceVec)).squaredNorm();
  actualEnergyDecrease =
      -f.energy.aggregationEnergy + previousE.aggregationEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "aggregation force convergence rate" << std::endl;

  // ==========================================================
  // ================   Self-avoidance       ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.selfAvoidanceForceVec));
  f.updateConfigurations(false);
  f.computeSelfAvoidanceEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.selfAvoidanceForceVec))
              .squaredNorm();
  actualEnergyDecrease =
      -f.energy.selfAvoidancePenalty + previousE.selfAvoidancePenalty;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.selfAvoidancePenalty <= previousE.selfAvoidancePenalty);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "self avoidance penalty = (expected - actual) / expected: "
      << difference_h / (abs(actualEnergyDecrease));

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h *
          f.forces.maskForce(toMatrix(f.forces.selfAvoidanceForceVec));
  f.updateConfigurations(false);
  f.computeSelfAvoidanceEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.selfAvoidanceForceVec))
          .squaredNorm();
  actualEnergyDecrease =
      -f.energy.selfAvoidancePenalty + previousE.selfAvoidancePenalty;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "self-avoidance force convergence rate" << std::endl;

  // ==========================================================
  // ================      Dirichlet         ==================
  // ==========================================================
  // test energy decrease and closeness
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      h * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec));
  f.updateConfigurations(false);
  f.computeDirichletEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec))
              .squaredNorm();
  actualEnergyDecrease = -f.energy.dirichletEnergy + previousE.dirichletEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.dirichletEnergy <= previousE.dirichletEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "line tension force = (expected - actual) / expected: "
      << difference_h / (abs(actualEnergyDecrease));

  // test convergence rate
  f.proteinDensity.raw() = current_proteinDensity;
  toMatrix(f.vpg->inputVertexPositions) =
      current_pos +
      stepFold * h *
          f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec));
  f.updateConfigurations(false);
  f.computeDirichletEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec))
          .squaredNorm();
  actualEnergyDecrease = -f.energy.dirichletEnergy + previousE.dirichletEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "dirichlet force convergence rate" << std::endl;
#pragma endregion force
#pragma region potential
  // ==========================================================
  // ================        Bending         ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.bendingPotential.raw());
  f.updateConfigurations(false);
  f.computeBendingEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.bendingPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.bendingEnergy + previousE.bendingEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.bendingEnergy <= previousE.bendingEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "bending potential: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h * f.forces.maskProtein(f.forces.bendingPotential.raw());
  f.updateConfigurations(false);
  f.computeBendingEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskProtein(f.forces.bendingPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.bendingEnergy + previousE.bendingEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "bending potential convergence rate" << std::endl;

  // ==========================================================
  // ================      Deviatoric        ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.deviatoricPotential.raw());
  f.updateConfigurations(false);
  f.computeDeviatoricEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskProtein(f.forces.deviatoricPotential.raw()).squaredNorm();
  actualEnergyDecrease =
      -f.energy.deviatoricEnergy + previousE.deviatoricEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.deviatoricEnergy <= previousE.deviatoricEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "deviatoric potential: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // no need to test convergence rate, expect exact since it is linear function

  // ==========================================================
  // ================      Adsorption        ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.adsorptionPotential.raw());
  f.updateConfigurations(false);
  f.computeAdsorptionEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskProtein(f.forces.adsorptionPotential.raw()).squaredNorm();
  actualEnergyDecrease =
      -f.energy.adsorptionEnergy + previousE.adsorptionEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.adsorptionEnergy <= previousE.adsorptionEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "adsorption force: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // no need to test convergence rate, expect exact since it is linear function

  // ==========================================================
  // ================      Aggregation       ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.aggregationPotential.raw());
  f.updateConfigurations(false);
  f.computeAggregationEnergy();
  expectedEnergyDecrease =
      h *
      f.forces.maskProtein(f.forces.aggregationPotential.raw()).squaredNorm();
  actualEnergyDecrease =
      -f.energy.aggregationEnergy + previousE.aggregationEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.aggregationEnergy <= previousE.aggregationEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "aggregation potential: expected vs. actual: "
      << expectedEnergyDecrease << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h * f.forces.maskProtein(f.forces.aggregationPotential.raw());
  f.updateConfigurations(false);
  f.computeAggregationEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskProtein(f.forces.aggregationPotential.raw()).squaredNorm();
  actualEnergyDecrease =
      -f.energy.aggregationEnergy + previousE.aggregationEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "aggregation potential convergence rate" << std::endl;

  // ==========================================================
  // ================      Dirichlet         ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.diffusionPotential.raw());
  f.updateConfigurations(false);
  f.computeDirichletEnergy();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.diffusionPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.dirichletEnergy + previousE.dirichletEnergy;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.dirichletEnergy <= previousE.dirichletEnergy);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "diffusion potential: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h * f.forces.maskProtein(f.forces.diffusionPotential.raw());
  f.updateConfigurations(false);
  f.computeDirichletEnergy();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskProtein(f.forces.diffusionPotential.raw()).squaredNorm();
  actualEnergyDecrease = -f.energy.dirichletEnergy + previousE.dirichletEnergy;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "diffusion potential convergence rate" << std::endl;

  // ==========================================================
  // ================     Interior Penalty   ==================
  // ==========================================================
  // test energy decrease and closeness
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      h * f.forces.maskProtein(f.forces.interiorPenaltyPotential.raw());
  f.updateConfigurations(false);
  f.computeProteinInteriorPenalty();
  expectedEnergyDecrease =
      h * f.forces.maskProtein(f.forces.interiorPenaltyPotential.raw())
              .squaredNorm();
  actualEnergyDecrease =
      -f.energy.proteinInteriorPenalty + previousE.proteinInteriorPenalty;
  difference_h = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_TRUE(f.energy.proteinInteriorPenalty <=
              previousE.proteinInteriorPenalty);
  EXPECT_TRUE(difference_h < tolerance * (abs(actualEnergyDecrease)))
      << "interior penalty: expected vs. actual: " << expectedEnergyDecrease
      << ", " << actualEnergyDecrease << std::endl;

  // test convergence rate
  toMatrix(f.vpg->inputVertexPositions) = current_pos;
  f.proteinDensity.raw() =
      current_proteinDensity +
      stepFold * h *
          f.forces.maskProtein(f.forces.interiorPenaltyPotential.raw());
  f.updateConfigurations(false);
  f.computeProteinInteriorPenalty();
  expectedEnergyDecrease =
      stepFold * h *
      f.forces.maskProtein(f.forces.interiorPenaltyPotential.raw())
          .squaredNorm();
  actualEnergyDecrease =
      -f.energy.proteinInteriorPenalty + previousE.proteinInteriorPenalty;
  difference_xh = abs(expectedEnergyDecrease - actualEnergyDecrease);
  EXPECT_NEAR(difference_xh / difference_h, pow(stepFold, expectRate),
              tolerance)
      << "interior penalty convergence rate" << std::endl;

#pragma endregion potential
};
} // namespace solver
} // namespace mem3dg
