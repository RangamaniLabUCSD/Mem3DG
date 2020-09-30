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

#include "mem3dg/solver/force.h"
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <cmath>

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getExternalForces() {

  auto externalPressure_e = gc::EigenMap<double, 3>(externalPressure);

  if (P.Kf != 0) {

    // a. FIND OUT THE CURRENT EXTERNAL PRESSURE MAGNITUDE BASED ON CURRENT GEOMETRY

    //auto &dist_e = heatMethodDistance(vpg, mesh.vertex(P.ptInd)).raw();
    //double stdDev = dist_e.maxCoeff() / P.conc;
    //externalPressureMagnitude =
    //    P.Kf / (stdDev * pow(M_PI * 2, 0.5)) *
    //    (-dist_e.array() * dist_e.array() / (2 * stdDev * stdDev)).exp();

    // b. APPLY EXTERNAL PRESSURE NORMAL TO THE SURFACE

    //auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
    //externalPressure_e = externalPressureMagnitude * vertexAngleNormal_e.row(P.ptInd);
    
    // c. ALTERNATIVELY, PRESSURE BASED ON INITIAL GEOMETRY + ALONG A FIXED DIRECTION, E.G. NEGATIVE Z DIRECTION

    Eigen::Matrix<double, 1, 3> zDir;
    zDir << 0.0, 0.0, -1.0;
    externalPressure_e =
        - externalPressureMagnitude * zDir * (vpg.inputVertexPositions[mesh.vertex(P.ptInd)].z - P.height);

  } else {
    externalPressure_e.setZero();
  }

}
} // namespace ddgsolver
