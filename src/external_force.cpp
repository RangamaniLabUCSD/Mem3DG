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
//     Padmini Rangmani (prangamani@eng.ucsd.edu)
//

#include <iostream>

#include <Eigen/Core>

#include "ddgsolver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getExternalForces() {
  auto externalForces_e = EigenMap<double, 3>(externalForces);
  if (P.extF != 0) {
    // std::cout << "external force::" << appliedForceMagnitude << std::endl;
    // externalForce_e = vertexAreaGradientNormal.array().colwise() *
    //	appliedForceMagnitude.array();
    auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);
    externalForces_e = appliedForceMagnitude * vertexAngleNormal_e.row(P.ptInd);
  } else {
    externalForces_e.setZero();
  }
}
} // namespace ddgsolver
