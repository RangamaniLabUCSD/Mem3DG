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

#include "mem3dg/mem3dg"
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
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> topologyMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vertexMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> refVertexMatrix;
  EigenVectorX1d proteinDensity;
  EigenVectorX3dr velocity;
  Parameters p;
  double h = 0.1;

  ForceTest() {
    // Create mesh and geometry objects
    std::tie(topologyMatrix, vertexMatrix) =
        getCylinderMatrix(1, 10, 10, 5, 0.3);
    std::tie(topologyMatrix, refVertexMatrix) = getCylinderMatrix(1, 10, 10);
    proteinDensity = Eigen::MatrixXd::Constant(vertexMatrix.rows(), 1, 1);
    velocity = Eigen::MatrixXd::Constant(vertexMatrix.rows(), 3, 0);

    p.variation.isShapeVariation = true;
    p.variation.isProteinVariation = true;
    p.variation.geodesicMask = -1;

    p.point.index = mem3dg::getVertexClosestToEmbeddedCoordinate(
        vertexMatrix, std::array<double, 3>{0, 0, 1});

    auto geodesicProteinDensity = [](double time, EigenVectorX1d meanCurvature,
                                     EigenVectorX1d geodesicDistance) {
      EigenVectorX1d proteinDensity;
      proteinDensity.resize(geodesicDistance.rows(), 1);
      for (std::size_t i = 0; i < geodesicDistance.rows(); i++) {
        proteinDensity[i] = tanhDistribution(geodesicDistance[i], 3, 1);
      }
      // tanhDistribution(proteinDensity, geodesicDistance, 10, 1);
      proteinDensity.array() *= 0.5;
      proteinDensity.array() += 0.2;
      return proteinDensity;
    };
    p.protein.form = geodesicProteinDensity;

    p.bending.alpha = 1;
    p.bending.dA0 = 4;
    p.bending.D = 0.01;
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

    p.aggregation.chi = -1e-4;

    p.entropy.xi = -1e-4;

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

    p.spring.Kst = 1e-3;
    p.spring.Kse = 1e-3;
    p.spring.Ksl = 1e-3;
  }
};

/**
 * @brief Test whether force is conservative: result need to be the same
 * when computed twice
 *
 */
TEST_F(ForceTest, ConservativeForcesTest) {
  // Instantiate system object
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, refVertexMatrix,
                           proteinDensity, velocity, p, 0);
  f.initialize(0, true);
  f.computeGeodesicDistance();
  f.prescribeProteinDensityDistribution();
  f.updateConfigurations();
  // First time calculation of force
  f.computeConservativeForcing();
  EigenVectorX3dr conservativeForceVec1 =
      toMatrix(f.forces.conservativeForceVec);
  EigenVectorX1d chemicalPotential1 = f.forces.chemicalPotential.raw();
  // Second time calculation of force
  f.computeConservativeForcing();
  EigenVectorX3dr conservativeForceVec2 =
      toMatrix(f.forces.conservativeForceVec);
  EigenVectorX1d chemicalPotential2 = f.forces.chemicalPotential.raw();

  // Comparison of 2 force calculations
  EXPECT_TRUE(conservativeForceVec1.isApprox(conservativeForceVec2));
  EXPECT_TRUE(chemicalPotential1.isApprox(chemicalPotential2));
  //   EXPECT_TRUE((conservativeForceVec1 - conservativeForceVec2).norm() <
  //   1e-12); EXPECT_TRUE((chemicalPotential1 - chemicalPotential2).norm() <
  //   1e-12); EXPECT_TRUE((regularizationForce1 - regularizationForce2).norm()
  //   < 1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to
 * 1. decrease in energy
 * 2. decrease in second order (or exact)
 */
TEST_F(ForceTest, ConsistentForceEnergy) {

  // initialize the system
  mem3dg::solver::System f(topologyMatrix, vertexMatrix, refVertexMatrix,
                           proteinDensity, velocity, p, 0);
  f.initialize(0, true);
  f.computeGeodesicDistance();
  f.prescribeProteinDensityDistribution();
  f.updateConfigurations();
  EXPECT_TRUE(f.testConservativeForcing(h));
};
} // namespace solver
} // namespace mem3dg
