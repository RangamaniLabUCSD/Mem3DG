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

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "mem3dg/solver/system.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::getBendingEnergy() {
  Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H - H0;
  E.BE = P.Kb * H_difference.transpose() * M * H_difference;
}

void System::getSurfaceEnergy() {
  double A_difference = surfaceArea - targetSurfaceArea;
  if (mesh.hasBoundary()) {
    E.sE = P.Ksg * A_difference;
  } else {
    E.sE = P.Ksg * A_difference * A_difference / targetSurfaceArea / 2 +
           P.lambdaSG * A_difference;
  }
}

void System::getPressureEnergy() {
  double V_difference = volume - refVolume * P.Vt;
  if (mesh.hasBoundary()) {
    E.pE = -P.Kv * V_difference;
  } else if (isReducedVolume) {
    E.pE = P.Kv * V_difference * V_difference / (refVolume * P.Vt) / 2 +
           P.lambdaV * V_difference;
  } else {
    E.pE = P.Kv * abs(log(P.Kv / volume / P.Pam));
  }
}

void System::getChemicalEnergy() {
  E.cE = (M * P.epsilon * proteinDensity.raw()).sum();
}

void System::getLineTensionEnergy() {
  E.lE = (P.eta * interArea * P.sharpness);
}

void System::getExternalForceEnergy() {
  E.exE = -rowwiseDotProduct(gc::EigenMap<double, 3>(externalPressure),
                             gc::EigenMap<double, 3>(vpg.inputVertexPositions))
               .sum();
}

void System::getKineticEnergy() {
  auto velocity = gc::EigenMap<double, 3>(vel);
  E.kE = 0.5 * (M * (velocity.array() * velocity.array()).matrix()).sum();
}

void System::getPotentialEnergy() {
  E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
}

void System::getFreeEnergy() {
  if (P.Kb != 0) {
    getBendingEnergy();
  }
  if (P.Ksg != 0) {
    getSurfaceEnergy();
  }
  if (P.Kv != 0) {
    getPressureEnergy();
  }
  if (isProtein) {
    getChemicalEnergy();
  }
  if (P.eta != 0) {
    getLineTensionEnergy();
  }
  if (P.Kf != 0) {
    getExternalForceEnergy();
  }
  getKineticEnergy();
  getPotentialEnergy();
  E.totalE = E.kE + E.potE;
}

void System::getL2ErrorNorm(
    Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure) {
  L2ErrorNorm =
      sqrt((M * rowwiseDotProduct(physicalPressure, physicalPressure)).sum());
}

void System::update_Vertex_positions() {
  vpg.refreshQuantities();

  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto positions = gc::EigenMap<double, 3>(vpg.inputVertexPositions);

  // initialize/update (inverse) mass and Laplacian matrix
  if (isTuftedLaplacian) {
    getTuftedLaplacianAndMass(M, L, mesh, vpg);
  } else {
    M = vpg.vertexLumpedMassMatrix;
    L = vpg.cotanLaplacian;
  }
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();

  // initialize/update distance from the point specified
  geodesicDistanceFromPtInd = heatSolver.computeDistance(mesh.vertex(ptInd));

  // initialize/update spontaneous curvature
  if (isProtein) {
    // proteinDensity.raw().setZero();
    // H0.setZero(mesh.nVertices(), 1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0 = (P.H0 * proteinDensitySq.array() / (1 + proteinDensitySq.array()))
             .matrix();
  } else if (P.H0 != 0) {
    tanhDistribution(vpg, H0, geodesicDistanceFromPtInd.raw(), P.sharpness,
                     P.r_H0);
    H0 *= P.H0;
    if (((H0.array() - (H0.sum() / mesh.nVertices())).matrix().norm() <
         1e-12)) {
      assert(P.eta == 0);
    }
  } else {
    H0.setZero(mesh.nVertices(), 1);
    assert(P.eta == 0);
  }

  // initialize/update mean curvature
  Eigen::Matrix<double, Eigen::Dynamic, 1> H_integrated =
      rowwiseDotProduct(L * positions / 2.0, vertexAngleNormal_e);
  H = M_inv * H_integrated;

  /// initialize/udate excess pressure
  volume = 0;
  for (gcs::Face f : mesh.faces()) {
    volume += signedVolumeFromFace(
        f, vpg, refVpg.inputVertexPositions[mesh.vertex(ptInd)]);
  }

  // initialize/update total surface area
  surfaceArea = vpg.faceAreas.raw().sum();

  // initialize/update intersection area
  interArea = 0.0;
  for (gcs::Vertex v : mesh.vertices()) {
    if ((H0[v.getIndex()] > (0.1 * P.H0)) &&
        (H0[v.getIndex()] < (0.9 * P.H0)) && (H[v.getIndex()] != 0)) {
      interArea += vpg.vertexDualAreas[v];
    }
  }

  // initialize/update external force
  getExternalPressure();

  // initialize/update the vertex position of the last iteration
  pastPositions = vpg.inputVertexPositions;

  // zero all forces
  gc::EigenMap<double, 3>(bendingPressure).setZero();
  gc::EigenMap<double, 3>(insidePressure).setZero();
  gc::EigenMap<double, 3>(capillaryPressure).setZero();
  gc::EigenMap<double, 3>(lineTensionPressure).setZero();
  gc::EigenMap<double, 3>(externalPressure).setZero();
  gc::EigenMap<double, 3>(regularizationForce).setZero();
  gc::EigenMap<double, 3>(dampingForce).setZero();
  gc::EigenMap<double, 3>(stochasticForce).setZero();
  chemicalPotential.raw().setZero();

  // zero all energy
  E = {0, 0, 0, 0, 0, 0, 0, 0, 0};
}
} // namespace mem3dg
