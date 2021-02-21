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

void System::computeBendingEnergy() {
  Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H.raw() - H0.raw();
  E.BE = P.Kb * H_difference.transpose() * M * H_difference;

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeSurfaceEnergy() {
  double A_difference = surfaceArea - targetSurfaceArea;
  if (mesh->hasBoundary()) {
    E.sE = P.Ksg * A_difference;
  } else {
    E.sE = P.Ksg * A_difference * A_difference / targetSurfaceArea / 2 +
           P.lambdaSG * A_difference;
  }
}

void System::computePressureEnergy() {
  double V_difference = volume - refVolume * P.Vt;
  if (mesh->hasBoundary()) {
    E.pE = -P.Kv * V_difference;
  } else if (isReducedVolume) {
    E.pE = P.Kv * V_difference * V_difference / (refVolume * P.Vt) / 2 +
           P.lambdaV * V_difference;
  } else {
    E.pE = P.Kv * (P.cam * volume - log(P.cam * volume) - 1);
  }
}

void System::computeChemicalEnergy() {
  E.cE = (M * P.epsilon * proteinDensity.raw()).sum();
}

void System::computeLineTensionEnergy() {
  E.lE = P.eta * (vpg->hodge1 * vpg->d0 * H0.raw().cwiseAbs()).sum();
}

void System::computeExternalForceEnergy() {
  E.exE = -rowwiseDotProduct(gc::EigenMap<double, 3>(externalPressure),
                             gc::EigenMap<double, 3>(vpg->inputVertexPositions))
               .sum();
}

void System::computeKineticEnergy() {
  // auto velocity = gc::EigenMap<double, 3>(vel);
  auto velocity =
      rowwiseDotProduct(gc::EigenMap<double, 3>(vel),
                        gc::EigenMap<double, 3>(vpg->inputVertexPositions));
  E.kE = 0.5 * (M * (velocity.array() * velocity.array()).matrix()).sum();
}

void System::computePotentialEnergy() {
  if (P.Kb != 0) {
    computeBendingEnergy();
  }
  if (P.Ksg != 0) {
    computeSurfaceEnergy();
  }
  if (P.Kv != 0) {
    computePressureEnergy();
  }
  if (isProtein) {
    computeChemicalEnergy();
  }
  if (P.eta != 0) {
    computeLineTensionEnergy();
  }
  if (P.Kf != 0) {
    computeExternalForceEnergy();
  }
  E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
}

void System::computeFreeEnergy() {
  // zero all energy
  E = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  computeKineticEnergy();
  computePotentialEnergy();
  E.totalE = E.kE + E.potE;
}

double
System::computeL2Norm(Eigen::Matrix<double, Eigen::Dynamic, 3> force) const {
  return sqrt((M * rowwiseDotProduct(force, force)).sum() / surfaceArea);

  // return sqrt(rowwiseDotProduct(force, force).sum()) / surfaceArea;

  // auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  // return (M * rowwiseDotProduct(pressure,
  // vertexAngleNormal_e).cwiseAbs()).sum() / surfaceArea;
}
} // namespace mem3dg
