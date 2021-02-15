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
  double h = 0.0002;
  bool isProtein = false, isLocalCurvature = false, isVertexShift = false,
       isReducedVolume = true;

  ForceCalculationTest() {
    // physical parameters
    p.H0 = 10;
    p.sharpness = 10;
    p.r_H0 = {0.5, 0.5};

    p.Vt = 0.6;
    p.cam = 0;

    p.pt = {1, 1, 1};
    p.Kf = 0;
    p.conc = 25;
    p.height = 0;

    p.Kb = 8.22e-5;
    p.Ksg = 0.01;
    p.Kv = 2e-2;
    p.eta = 0;

    p.epsilon = 15e-5;
    p.Bc = 40;

    p.Kse = 0;
    p.Ksl = 0;
    p.Kst = 6;

    p.gamma = 0;
    p.temp = 0;
    p.sigma = 0;

    p.radius = 3;

    // Create mesh and geometry objects
    std::tie(ptrMesh, ptrVpg) = icosphere(2, 1);
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
                   isReducedVolume, isProtein, isLocalCurvature, isVertexShift);

  // First time calculation of force
  f.getAllForces();
  EigenVectorX3D bendingPressure1 = gc::EigenMap<double, 3>(f.bendingPressure),
                 insidePressure1 = f.insidePressure * gc::EigenMap<double, 3>(
                                                          f.vpg->vertexNormals),
                 capillaryPressure1 =
                     gc::EigenMap<double, 3>(f.capillaryPressure),
                 regularizationForce1 =
                     gc::EigenMap<double, 3>(f.regularizationForce),
                 lineTensionPressure1 =
                     gc::EigenMap<double, 3>(f.lineTensionPressure),
                 externalPressure1 =
                     gc::EigenMap<double, 3>(f.externalPressure);
  EigenVectorX1D chemicalPotential1 = f.chemicalPotential.raw();

  // Second time calculation of force
  f.getAllForces();
  EigenVectorX3D bendingPressure2 = gc::EigenMap<double, 3>(f.bendingPressure),
                 insidePressure2 = f.insidePressure * gc::EigenMap<double, 3>(
                                                          f.vpg->vertexNormals),
                 capillaryPressure2 =
                     gc::EigenMap<double, 3>(f.capillaryPressure),
                 regularizationForce2 =
                     gc::EigenMap<double, 3>(f.regularizationForce),
                 lineTensionPressure2 =
                     gc::EigenMap<double, 3>(f.lineTensionPressure),
                 externalPressure2 =
                     gc::EigenMap<double, 3>(f.externalPressure);
  EigenVectorX1D chemicalPotential2 = f.chemicalPotential.raw();

  // Comparison of 2 force calculations
  ASSERT_TRUE((bendingPressure1 - bendingPressure2).norm() < 1e-12);
  ASSERT_TRUE((capillaryPressure1 - capillaryPressure2).norm() < 1e-12);
  ASSERT_TRUE((insidePressure1 - insidePressure2).norm() < 1e-12);
  ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
  ASSERT_TRUE((externalPressure1 - externalPressure2).norm() < 1e-12);
  ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  ASSERT_TRUE((lineTensionPressure1 - lineTensionPressure2).norm() < 1e-12);
};

/**
 * @brief Test whether one-pass force computation is consistent with the
 * individual component force computation
 *
 */
TEST_F(ForceCalculationTest, OnePassVsReferenceForce) {
  // Instantiate system object
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   isReducedVolume, isProtein, isLocalCurvature, isVertexShift);

  // Get forces in one-pass
  f.getAllForces();
  EigenVectorX3D bendingPressure1 = gc::EigenMap<double, 3>(f.bendingPressure),
                 insidePressure1 = f.insidePressure * gc::EigenMap<double, 3>(
                                                          f.vpg->vertexNormals),
                 capillaryPressure1 =
                     gc::EigenMap<double, 3>(f.capillaryPressure),
                 regularizationForce1 =
                     gc::EigenMap<double, 3>(f.regularizationForce),
                 lineTensionPressure1 =
                     gc::EigenMap<double, 3>(f.lineTensionPressure),
                 externalPressure1 =
                     gc::EigenMap<double, 3>(f.externalPressure);
  // chemicalPotential1 = f.chemicalPotential.raw();

  // Get force individually
  f.getBendingPressure();
  f.getCapillaryPressure();
  f.getInsidePressure();
  f.getRegularizationForce();
  f.getLineTensionPressure();
  f.getExternalPressure();
  //   f.getChemicalPotential();
  EigenVectorX3D bendingPressure2 = gc::EigenMap<double, 3>(f.bendingPressure),
                 insidePressure2 = f.insidePressure * gc::EigenMap<double, 3>(
                                                          f.vpg->vertexNormals),
                 capillaryPressure2 =
                     gc::EigenMap<double, 3>(f.capillaryPressure),
                 regularizationForce2 =
                     gc::EigenMap<double, 3>(f.regularizationForce),
                 lineTensionPressure2 =
                     gc::EigenMap<double, 3>(f.lineTensionPressure),
                 externalPressure2 =
                     gc::EigenMap<double, 3>(f.externalPressure);
  //   EigenVectorX1D chemicalPotential2 = f.chemicalPotential.raw();

  // Comparison of two force calculations
  ASSERT_TRUE((bendingPressure1 - bendingPressure2).norm() < 1e-12);
  ASSERT_TRUE((capillaryPressure1 - capillaryPressure2).norm() < 1e-12);
  ASSERT_TRUE((insidePressure1 - insidePressure2).norm() < 1e-12);
  ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
  ASSERT_TRUE((externalPressure1 - externalPressure2).norm() < 1e-12);
  //   ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  ASSERT_TRUE((lineTensionPressure1 - lineTensionPressure2).norm() < 1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to decrease in
 * energy
 *
 */
TEST_F(ForceCalculationTest, ConsistentForceEnergy) {
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   isReducedVolume, isProtein, isLocalCurvature, isVertexShift);
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  //   Energy E_pre{f.E.totalE, f.E.kE, f.E.potE, f.E.BE, f.E.sE,
  //                f.E.pE,     f.E.cE, f.E.lE,   f.E.exE};
  f.getFreeEnergy();
  Energy E_pre{f.E};
  Energy E_aft;
  //   double totalE_pre = f.E.totalE, kE_pre = f.E.kE, potE_pre = f.E.potE,
  //          BE_pre = f.E.BE, sE_pre = f.E.sE, pE_pre = f.E.pE, cE_pre =
  //          f.E.cE, lE_pre = f.E.lE, exE = f.E.exE;
  for (size_t i = 0; i < 50; i++) {
    f.getBendingPressure();
    vel_e = rowwiseScaling(f.mask.cast<double>(),
                           gc::EigenMap<double, 3>(f.bendingPressure));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.getFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.BE <= E_pre.BE);
    E_pre = E_aft;

    f.getCapillaryPressure();
    vel_e = rowwiseScaling(f.mask.cast<double>(),
                           gc::EigenMap<double, 3>(f.capillaryPressure));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.getFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.sE <= E_pre.sE);
    E_pre = E_aft;

    f.getInsidePressure();
    vel_e = rowwiseScaling(f.mask.cast<double>(),
                           f.insidePressure *
                               gc::EigenMap<double, 3>(f.vpg->vertexNormals));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.getFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.pE <= E_pre.pE);
    E_pre = E_aft;

    f.getExternalPressure();
    vel_e = rowwiseScaling(f.mask.cast<double>(),
                           gc::EigenMap<double, 3>(f.externalPressure));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.getFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.exE <= E_pre.exE);
    E_pre = E_aft;

    f.getRegularizationForce();
    vel_e = rowwiseScaling(f.mask.cast<double>(),
                           gc::EigenMap<double, 3>(f.externalPressure));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.getFreeEnergy();
    E_aft = f.E;
    E_pre = E_aft;
  }

}; // namespace ddgsolver
} // namespace mem3dg
