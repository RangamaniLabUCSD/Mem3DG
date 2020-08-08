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

#include <Eigen/Core>

#include "ddgsolver/force.h"

#include <cmath>

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getExternalForces() {

  auto externalForces_e = EigenMap<double, 3>(externalForces);

  if (P.Kf != 0) {

    //// Find out the applied force magnitude based on current geometry
    //auto &dist_e = heatMethodDistance(vpg, mesh.vertex(P.ptInd)).raw();
    //double stdDev = dist_e.maxCoeff() / P.conc;
    //appliedForceMagnitude =
    //    P.Kf / (stdDev * pow(M_PI * 2, 0.5)) *
    //    (-dist_e.array() * dist_e.array() / (2 * stdDev * stdDev)).exp();

    //// Apply external forces always normal to the surface 
    //auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);
    //externalForces_e = appliedForceMagnitude * vertexAngleNormal_e.row(P.ptInd);
    
    // Alternatively, applied force maganitude based on initial geometry along a fixed direction, e.g. negative z direction
    Eigen::Matrix<double, 1, 3> zDir;
    zDir << 0.0, 0.0, -1.0;
    externalForces_e =
        - appliedForceMagnitude * zDir * (vpg.inputVertexPositions[mesh.vertex(P.ptInd)].z - P.height);
    ;

  } else {
    externalForces_e.setZero();
  }

}
} // namespace ddgsolver
