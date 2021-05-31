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
  double h = 0.0002;

  ForceCalculationTest() {
    // physical parameters
    p.H0c = 10;
    p.protein0.resize(1, 1);
    p.protein0 << 0.5;

    p.Vt = 0.6;
    p.cam = 0;

    p.pt.resize(2, 1);
    p.pt << 1, 1;
    p.Kf = 0;
    p.conc = 25;
    p.height = 0;

    p.Kb = 8.22e-5;
    p.Ksg = 0.01;
    p.A_res = 0;
    p.Kv = 2e-2;
    p.V_res = 0;
    p.eta = 0;

    p.epsilon = 15e-5;
    p.Bc = 40;

    p.Kse = 0;
    p.Ksl = 0;
    p.Kst = 6;

    p.gamma = 0;
    p.temp = 0;

    p.radius = 3;

    o.isVertexShift = false;
    o.isProteinVariation = false;
    o.isReducedVolume = true;
    o.isEdgeFlip = false;
    o.isSplitEdge = false;
    o.isCollapseEdge = false;

    // Create mesh and geometry objects
    std::tie(ptrMesh, ptrVpg) = icosphere(1, 3);
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

  // First time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX1D bendingForce1 = f.F.bendingForce.raw(),
                 osmoticForce1 = f.F.osmoticForce.raw(),
                 capillaryForce1 = f.F.capillaryForce.raw(),
                 lineTensionForce1 = f.F.lineCapillaryForce.raw(),
                 externalForce1 = f.F.externalForce.raw(),
                 chemicalPotential1 = f.F.chemicalPotential.raw();
  EigenVectorX3D regularizationForce1 =
      gc::EigenMap<double, 3>(f.F.regularizationForce);

  // Second time calculation of force
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX1D bendingForce2 = f.F.bendingForce.raw(),
                 osmoticForce2 = f.F.osmoticForce.raw(),
                 capillaryForce2 = f.F.capillaryForce.raw(),
                 lineTensionForce2 = f.F.lineCapillaryForce.raw(),
                 externalForce2 = f.F.externalForce.raw(),
                 chemicalPotential2 = f.F.chemicalPotential.raw();
  EigenVectorX3D regularizationForce2 =
      gc::EigenMap<double, 3>(f.F.regularizationForce);

  // Comparison of 2 force calculations
  ASSERT_TRUE((bendingForce1 - bendingForce2).norm() < 1e-12);
  ASSERT_TRUE((capillaryForce1 - capillaryForce2).norm() < 1e-12);
  ASSERT_TRUE((osmoticForce1 - osmoticForce2).norm() < 1e-12);
  ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
  ASSERT_TRUE((externalForce1 - externalForce2).norm() < 1e-12);
  ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  ASSERT_TRUE((lineTensionForce1 - lineTensionForce2).norm() < 1e-12);
};

/**
 * @brief Test whether one-pass force computation is consistent with the
 * individual component force computation
 *
 */
TEST_F(ForceCalculationTest, OnePassVsReferenceForce) {
  // Instantiate system object
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   o);

  // Get forces in one-pass
  f.computePhysicalForces();
  f.computeRegularizationForce();
  EigenVectorX1D bendingForce1 = f.F.bendingForce.raw(),
                 osmoticForce1 = f.F.osmoticForce.raw(),
                 capillaryForce1 = f.F.capillaryForce.raw(),
                 lineTensionForce1 = f.F.lineCapillaryForce.raw(),
                 externalForce1 = f.F.externalForce.raw();
  //  chemicalPotential1 = f.chemicalPotential.raw();
  EigenVectorX3D regularizationForce1 =
      gc::EigenMap<double, 3>(f.F.regularizationForce);

  // Get force individually
  f.computeBendingForce();
  f.computeCapillaryForce();
  f.computeOsmoticForce();
  f.computeRegularizationForce();
  f.computeLineCapillaryForce();
  f.computeExternalForce();
  //   f.computeChemicalPotential();
  EigenVectorX1D bendingForce2 = f.F.bendingForce.raw(),
                 osmoticForce2 = f.F.osmoticForce.raw(),
                 capillaryForce2 = f.F.capillaryForce.raw(),
                 lineTensionForce2 = f.F.lineCapillaryForce.raw(),
                 externalForce2 = f.F.externalForce.raw();
  //  chemicalPotential2 = f.chemicalPotential.raw();
  EigenVectorX3D regularizationForce2 =
      gc::EigenMap<double, 3>(f.F.regularizationForce);

  // Comparison of two force calculations
  ASSERT_TRUE((bendingForce1 - bendingForce2).norm() < 1e-12);
  ASSERT_TRUE((capillaryForce1 - capillaryForce2).norm() < 1e-12);
  ASSERT_TRUE((osmoticForce1 - osmoticForce2).norm() < 1e-12);
  ASSERT_TRUE((regularizationForce1 - regularizationForce2).norm() < 1e-12);
  ASSERT_TRUE((externalForce1 - externalForce2).norm() < 1e-12);
  //   ASSERT_TRUE((chemicalPotential1 - chemicalPotential2).norm() < 1e-12);
  ASSERT_TRUE((lineTensionForce1 - lineTensionForce2).norm() < 1e-12);
};

/**
 * @brief Test whether integrating with the force will lead to decrease in
 * energy
 *
 */
TEST_F(ForceCalculationTest, ConsistentForceEnergy) {
  mem3dg::System f(std::move(ptrMesh), std::move(ptrVpg), std::move(ptrVpg), p,
                   o);
  auto vel_e = gc::EigenMap<double, 3>(f.vel);
  auto pos_e = gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  //   Energy E_pre{f.E.totalE, f.E.kE, f.E.potE, f.E.BE, f.E.sE,
  //                f.E.pE,     f.E.cE, f.E.lE,   f.E.exE};
  f.computeFreeEnergy();
  Energy E_pre{f.E};
  Energy E_aft;
  //   double totalE_pre = f.E.totalE, kE_pre = f.E.kE, potE_pre = f.E.potE,
  //          BE_pre = f.E.BE, sE_pre = f.E.sE, pE_pre = f.E.pE, cE_pre =
  //          f.E.cE, lE_pre = f.E.lE, exE = f.E.exE;
  for (size_t i = 0; i < 50; i++) {
    f.computeBendingForce();
    vel_e = f.F.mask(rowwiseScaling(f.F.bendingForce.raw(),
                                    EigenMap<double, 3>(f.vpg->vertexNormals)));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.BE <= E_pre.BE);
    E_pre = E_aft;

    f.computeCapillaryForce();
    vel_e = f.F.mask(rowwiseScaling(f.F.capillaryForce.raw(),
                                    EigenMap<double, 3>(f.vpg->vertexNormals)));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.sE <= E_pre.sE);
    E_pre = E_aft;

    f.computeOsmoticForce();
    vel_e = f.F.mask(rowwiseScaling(f.F.osmoticForce.raw(),
                                    EigenMap<double, 3>(f.vpg->vertexNormals)));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.pE <= E_pre.pE);
    E_pre = E_aft;

    f.computeExternalForce();
    vel_e = f.F.mask(rowwiseScaling(f.F.externalForce.raw(),
                                    EigenMap<double, 3>(f.vpg->vertexNormals)));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    E_aft = f.E;
    ASSERT_TRUE(E_aft.exE <= E_pre.exE);
    E_pre = E_aft;

    f.computeRegularizationForce();
    vel_e = f.F.mask(gc::EigenMap<double, 3>(f.F.regularizationForce));
    pos_e += vel_e * h;
    f.updateVertexPositions();
    f.computeFreeEnergy();
    E_aft = f.E;
    E_pre = E_aft;
  }

}; // namespace ddgsolver
} // namespace mem3dg
