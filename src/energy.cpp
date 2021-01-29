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

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // vpg.vertexGaussianCurvatures.raw()).sum();
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
    E.pE = P.Kv * (P.cam * volume - log(P.cam * volume) - 1);
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
  // auto velocity = gc::EigenMap<double, 3>(vel);
  auto velocity =
      rowwiseDotProduct(gc::EigenMap<double, 3>(vel),
                        gc::EigenMap<double, 3>(vpg.inputVertexPositions));
  E.kE = 0.5 * (M * (velocity.array() * velocity.array()).matrix()).sum();
}

void System::getPotentialEnergy() {
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
  E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
}

void System::getFreeEnergy() {
  // zero all energy
  E = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  getKineticEnergy();
  getPotentialEnergy();
  E.totalE = E.kE + E.potE;
}

void System::getL2ErrorNorm(Eigen::Matrix<double, Eigen::Dynamic, 3> pressure) {
  L2ErrorNorm =
      sqrt((M * rowwiseDotProduct(M * pressure, M * pressure)).sum() / surfaceArea);
}

double System::getL2Norm(gcs::VertexData<gc::Vector3> pressure_data) const {
  auto pressure = gc::EigenMap<double, 3>(pressure_data);
  return sqrt((M * rowwiseDotProduct(M * pressure, M * pressure)).sum() / surfaceArea);
}

} // namespace mem3dg
