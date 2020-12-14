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

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/integrator.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>

namespace ddgsolver {
namespace integration {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

double
getL2ErrorNorm(Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure) {
  return sqrt((rowwiseDotProduct(physicalPressure, physicalPressure)).sum());
}

std::tuple<double, double, double, double, double, double, double, double>
getFreeEnergy(Force &f) {
  double bE = 0.0, sE = 0.0, pE = 0.0, kE = 0.0, exE = 0.0, cE = 0.0, lE = 0.0,
         totalE = 0.0;

  if (f.mesh.hasBoundary()) {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = f.H - f.H0;
    double A_difference = f.surfaceArea - f.targetSurfaceArea;
    double V_difference = f.volume - f.refVolume * f.P.Vt;

    bE = (f.P.Kb * f.M *
          (f.mask.cast<double>().array() * H_difference.array() *
           H_difference.array())
              .matrix())
             .sum();
    sE = f.P.Ksg * A_difference;
    pE = -f.P.Kv * V_difference;

    auto vel = gc::EigenMap<double, 3>(f.vel);
    kE = 0.5 * (f.M * (vel.array() * vel.array()).matrix()).sum();

    if (f.isProtein) {
      cE = (f.M * f.P.epsilon * f.proteinDensity.raw()).sum();
    }

    lE = (f.P.eta * f.interArea * f.P.sharpness);

    exE =
        -rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                           gc::EigenMap<double, 3>(f.vpg.inputVertexPositions))
             .sum();

    totalE = bE + sE + pE + kE + cE + lE + exE;

  } else {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = f.H - f.H0;
    double A_difference = f.surfaceArea - f.targetSurfaceArea;
    double V_difference = f.volume - f.refVolume * f.P.Vt;

    bE = (f.P.Kb * f.M *
          (f.mask.cast<double>().array() * H_difference.array() *
           H_difference.array())
              .matrix())
             .sum();
    sE = f.P.Ksg * A_difference * A_difference / f.targetSurfaceArea / 2;
    pE = f.P.Kv * V_difference * V_difference / (f.refVolume * f.P.Vt) / 2;

    auto vel = gc::EigenMap<double, 3>(f.vel);
    kE = 0.5 * (f.M * (vel.array() * vel.array()).matrix()).sum();

    if (f.isProtein) {
      cE = (f.M * f.P.epsilon * f.proteinDensity.raw()).sum();
    }

    lE = (f.P.eta * f.interArea * f.P.sharpness);

    exE =
        -rowwiseDotProduct(gc::EigenMap<double, 3>(f.externalPressure),
                           gc::EigenMap<double, 3>(f.vpg.inputVertexPositions))
             .sum();

    totalE = bE + sE + pE + kE + cE + lE + exE;
  }

  std::tuple<double, double, double, double, double, double, double, double> output(
      totalE, bE, sE, pE, kE, cE, lE, exE);

  return output;
}
} // namespace integration
} // namespace ddgsolver
