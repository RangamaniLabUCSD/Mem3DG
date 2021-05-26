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
#include "mem3dg/solver/util.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

class ForceCalculationTest : public testing::Test {
protected:
  // initialize mesh and vpg
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  Parameters p;
  Options o;
  double h = 0.0001;

  ForceCalculationTest() {
    // options
    o.isProteinVariation = true;
    o.isReducedVolume = true;
    o.isConstantOsmoticPressure = false;
    o.isConstantSurfaceTension = false;
    o.isHeterogeneous = false;
    o.isEdgeFlip = true;
    o.isSplitEdge = true;
    o.isCollapseEdge = true;
    o.isVertexShift = true;
    o.isRefMesh = false;
    o.isFloatVertex = true;
    o.isLaplacianMeanCurvature = false;
    o.boundaryConditionType = "none";

    // physical parameters
    p.pt = {0, 0, 1};
    p.r_heter = {-1, -1};
    // bending //
    p.Kb = 8.22e-5;
    p.Kbc = 0; // 8.22e-4 DEFINITION OF LARGE AND SMALL VALUE;
    p.H0 = 10;
    // surface tension //
    p.Ksg = 1e-3 * 100;
    p.A_res = 0;
    p.epsilon = -1e-2;
    // osmotic force //
    p.Kv = 0.02;
    p.V_res = 4 / 3 * 3.14 * 1 * 1 * 1 * 0;
    p.Vt = 0.6;
    p.cam = -1;
    // regularization //
    p.Kst = 0; // 2e-6;
    p.Ksl = 0;
    p.Kse = 0;
    // protein binding //
    p.Bc = 0.5;
    // line tension //
    p.eta = 0.1;
    // DPD //
    p.gamma = 0;
    p.temp = 0;
    // Rarely used //
    p.Kf = 0;
    p.conc = -1;
    p.height = 0;
    p.radius = 100000;
    p.lambdaSG = 0;
    p.lambdaV = 0;

    // Create mesh and geometry objects
    std::tie(ptrMesh, ptrVpg) = icosphere(3, 1);
  }
};

/**
 * @brief Test whether passive force is conservative: result need to be the same
 * when computed twice
 *
 */
TEST_F(ForceCalculationTest, ConsistentForcesTest) {
  // Instantiate system object
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   o);

//   // First time calculation of force
//   f.computePhysicalForces();
//   f.computeRegularizationForce();
//   EigenVectorX3D vectorForces1 = f.F.toMatrix(f.F.vectorForces);
//   EigenVectorX1D chemicalPotential1 = f.F.toMatrix(f.F.chemicalPotential);
//   EigenVectorX3D regularizationForce1 = f.F.toMatrix(f.F.regularizationForce);

//   // Second time calculation of force
//   f.computePhysicalForces();
//   f.computeRegularizationForce();
//   EigenVectorX3D vectorForces2 = f.F.toMatrix(f.F.vectorForces);
//   EigenVectorX1D chemicalPotential2 = f.F.toMatrix(f.F.chemicalPotential);
//   EigenVectorX3D regularizationForce2 = f.F.toMatrix(f.F.regularizationForce);

//   // Comparison of 2 force calculations
//   ASSERT_TRUE((vectorForces1 - vectorForces2).norm() < 1e-12);
//   ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
//   ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to decrease in
 * energy
 *
 */
TEST_F(ForceCalculationTest, ConsistentForceEnergy) {
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   o);
//   auto vel_e = gc::EigenMap<double, 3>(f.vel);
//   auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
//   const EigenVectorX3D current_pos = f.F.toMatrix(f.vpg->inputVertexPositions);
//   const EigenVectorX1D current_proteinDensity = f.F.toMatrix(f.proteinDensity);
//   f.updateVertexPositions();
//   f.computeFreeEnergy();
//   Energy previousE{f.E};
//   f.computePhysicalForces();

//   // bending force
//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.toMatrix(f.F.bendingForceVec));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.BE <= previousE.BE);
//   std::cout << "With only bending force, actual energy decrease: "
//             << f.E.BE - previousE.BE << " vs. expected decrease"
//             << -h * (f.F.mask(f.F.toMatrix(f.F.bendingForceVec)).array() *
//                      f.F.mask(f.F.toMatrix(f.F.bendingForceVec)).array())
//                         .sum()
//             << std::endl;

//   f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
//   f.proteinDensity.raw() =
//       current_proteinDensity + h * f.P.Bc * f.F.bendingPotential.raw();
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.BE <= previousE.BE);
//   std::cout << "With only bending potential, actual energy decrease: "
//             << f.E.BE - previousE.BE << " from " << previousE.BE << " to "
//             << f.E.BE << ", expected dBE: "
//             << -h * f.P.Bc *
//                    (f.F.bendingPotential.raw().array() *
//                     f.F.bendingPotential.raw().array())
//                        .sum()
//             << std::endl;

//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.toMatrix(f.F.capillaryForceVec));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.sE <= previousE.sE);
//   std::cout << "With only capillary force, actual energy decrease: "
//             << f.E.sE - previousE.sE << " from " << previousE.sE << " to "
//             << f.E.sE << ", expected dsE: "
//             << -h * (f.F.mask(f.F.toMatrix(f.F.capillaryForceVec)).array() *
//                      f.F.mask(f.F.toMatrix(f.F.capillaryForceVec)).array())
//                         .sum()
//             << std::endl;

//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.toMatrix(f.F.osmoticForceVec));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.pE <= previousE.pE);
//   std::cout << "With only osmotic force, actual energy decrease: "
//             << f.E.pE - previousE.pE << " from " << previousE.pE << " to "
//             << f.E.pE << ", expected dpE: "
//             << -h * (f.F.mask(f.F.toMatrix(f.F.osmoticForceVec)).array() *
//                      f.F.mask(f.F.toMatrix(f.F.osmoticForceVec)).array())
//                         .sum()
//             << std::endl;

//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.aE <= previousE.aE);
//   std::cout << "With only adsorption force, actual energy decrease: "
//             << f.E.aE - previousE.aE << " from " << previousE.aE << " to "
//             << f.E.aE << ", expected daE: "
//             << -h * (f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec)).array() *
//                      f.F.mask(f.F.toMatrix(f.F.adsorptionForceVec)).array())
//                         .sum()
//             << std::endl;

//   f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
//   f.proteinDensity.raw() =
//       current_proteinDensity + h * f.P.Bc * f.F.adsorptionPotential.raw();
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.aE <= previousE.aE);
//   std::cout << "With only adsorption potential, actual energy decrease: "
//             << f.E.aE - previousE.aE << " from " << previousE.aE << " to "
//             << f.E.aE << ", expected dBE: "
//             << -h * f.P.Bc *
//                    (f.F.adsorptionPotential.raw().array() *
//                     f.F.adsorptionPotential.raw().array())
//                        .sum()
//             << std::endl;

//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.dE <= previousE.dE);
//   std::cout << "With only line tension force, actual energy decrease: "
//             << f.E.dE - previousE.dE << " from " << previousE.dE << " to "
//             << f.E.dE << ", expected ddE: "
//             << -h * (f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec)).array() *
//                      f.F.mask(f.F.toMatrix(f.F.lineCapillaryForceVec)).array())
//                         .sum()
//             << std::endl;

//   f.F.toMatrix(f.vpg->inputVertexPositions) = current_pos;
//   f.proteinDensity.raw() =
//       current_proteinDensity + h * f.P.Bc * f.F.diffusionPotential.raw();
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.dE <= previousE.dE);
//   std::cout << "With only diffusion potential, actual energy decrease: "
//             << f.E.dE - previousE.dE << " from " << previousE.dE << " to "
//             << f.E.dE << ", expected ddE: "
//             << -h * f.P.Bc *
//                    (f.F.diffusionPotential.raw().array() *
//                     f.F.diffusionPotential.raw().array())
//                        .sum()
//             << std::endl;

//   f.proteinDensity.raw() = current_proteinDensity;
//   f.F.toMatrix(f.vpg->inputVertexPositions) =
//       current_pos + h * f.F.mask(f.F.addNormal(f.F.externalForce.raw()));
//   f.updateVertexPositions(false);
//   f.computeFreeEnergy();
//   ASSERT_TRUE(f.E.exE <= previousE.exE);
//   std::cout << "With only external force, actual energy decrease: "
//             << f.E.exE - previousE.exE << " from " << previousE.exE << " to "
//             << f.E.exE << ", expected dexE: "
//             << -h * (f.F.mask(f.F.addNormal(f.F.externalForce.raw())).array() *
//                      f.F.mask(f.F.addNormal(f.F.externalForce.raw())).array())
//                         .sum()
//             << std::endl;

  // for (size_t i = 0; i < 50; i++) {
  //   f.computeBendingForce();
  //   vel_e = f.F.mask(rowwiseScaling(f.F.bendingForce.raw(),
  //                                   EigenMap<double,
  //                                   3>(f.vpg->vertexNormals)));
  //   pos_e += vel_e * h;
  //   f.updateVertexPositions();
  //   f.computeFreeEnergy();
  //   E_aft = f.E;
  //   ASSERT_TRUE(E_aft.BE <= previousE.BE);
  //   previousE = E_aft;

  //   f.computeCapillaryForce();
  //   vel_e = f.F.mask(rowwiseScaling(f.F.capillaryForce.raw(),
  //                                   EigenMap<double,
  //                                   3>(f.vpg->vertexNormals)));
  //   pos_e += vel_e * h;
  //   f.updateVertexPositions();
  //   f.computeFreeEnergy();
  //   E_aft = f.E;
  //   ASSERT_TRUE(E_aft.sE <= previousE.sE);
  //   previousE = E_aft;

  //   f.computeOsmoticForce();
  //   vel_e = f.F.mask(rowwiseScaling(f.F.osmoticForce.raw(),
  //                                   EigenMap<double,
  //                                   3>(f.vpg->vertexNormals)));
  //   pos_e += vel_e * h;
  //   f.updateVertexPositions();
  //   f.computeFreeEnergy();
  //   E_aft = f.E;
  //   ASSERT_TRUE(E_aft.pE <= previousE.pE);
  //   previousE = E_aft;

  //   f.computeExternalForce();
  //   vel_e = f.F.mask(rowwiseScaling(f.F.externalForce.raw(),
  //                                   EigenMap<double,
  //                                   3>(f.vpg->vertexNormals)));
  //   pos_e += vel_e * h;
  //   f.updateVertexPositions();
  //   f.computeFreeEnergy();
  //   E_aft = f.E;
  //   ASSERT_TRUE(E_aft.exE <= previousE.exE);
  //   previousE = E_aft;

  //   f.computeRegularizationForce();
  //   vel_e = f.F.mask(gc::EigenMap<double, 3>(f.F.regularizationForce));
  //   pos_e += vel_e * h;
  //   f.updateVertexPositions();
  //   f.computeFreeEnergy();
  //   E_aft = f.E;
  //   previousE = E_aft;
  // }

}; // namespace ddgsolver
} // namespace mem3dg
