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

#include "mem3dg/solver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::getL2ErrorNorm(
    Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure) {
  L2ErrorNorm =
      sqrt((M * rowwiseDotProduct(physicalPressure, physicalPressure)).sum());
}

void System::getFreeEnergy() {

  if (mesh.hasBoundary()) {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H - H0;
    double A_difference = surfaceArea - targetSurfaceArea;
    double V_difference = volume - refVolume * P.Vt;

    E.BE = P.Kb * H_difference.transpose() * M * H_difference;
    E.sE = P.Ksg * A_difference;
    E.pE = -P.Kv * V_difference;

    auto velocity = gc::EigenMap<double, 3>(vel);
    E.kE = 0.5 * (M * (velocity.array() * velocity.array()).matrix()).sum();

    if (isProtein) {
      E.cE = (M * P.epsilon * proteinDensity.raw()).sum();
    }

    E.lE = (P.eta * interArea * P.sharpness);

    E.exE =
        -rowwiseDotProduct(gc::EigenMap<double, 3>(externalPressure),
                           gc::EigenMap<double, 3>(vpg.inputVertexPositions))
             .sum();
    E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
    E.totalE = E.potE + E.kE;

  } else {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H - H0;
    double A_difference = surfaceArea - targetSurfaceArea;
    double V_difference = volume - refVolume * P.Vt;

    E.BE = P.Kb * H_difference.transpose() * M * H_difference;
    E.sE = P.Ksg * A_difference * A_difference / targetSurfaceArea / 2 +
           P.lambdaSG * A_difference;
    E.pE = P.Kv * V_difference * V_difference / (refVolume * P.Vt) / 2 +
           P.lambdaV * V_difference;

    auto velocity = gc::EigenMap<double, 3>(vel);
    E.kE = 0.5 * (M * (velocity.array() * velocity.array()).matrix()).sum();

    if (isProtein) {
      E.cE = (M * P.epsilon * proteinDensity.raw()).sum();
    }

    E.lE = (P.eta * interArea * P.sharpness);

    E.exE =
        -rowwiseDotProduct(gc::EigenMap<double, 3>(externalPressure),
                           gc::EigenMap<double, 3>(vpg.inputVertexPositions))
             .sum();

    E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
    E.totalE = E.potE + E.kE;
  }
}
} // namespace ddgsolver
