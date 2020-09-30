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
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>

#include <iostream>

namespace ddgsolver {
namespace integration {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

double getL2ErrorNorm(
    Eigen::SparseMatrix<double> M,
    Eigen::Matrix<double, Eigen::Dynamic, 3>
        physicalPressure){

  return sqrt(
      (M * rowwiseDotProduct(physicalPressure, physicalPressure)).sum());
}

std::tuple<double, double, double, double, double> getFreeEnergy(Force &f) {
  // comment: this may not be useful, the convergence can be be tested by
  // checking its derivative which is the forces excluding the DPD forces. The
  // energy trajectory of could actually numerically integrated by post
  // processing after saving all forces during the iterations.

  // auto positions = gc::EigenMap<double, 3>(vpg.inputVertexPositions);
  // auto A = f.L.transpose() * f.M_inv * f.L; //could further cached since the
  // same for all time pt for (size_t i = 0; i < positions.rows(); i++) {
  // for
  //(size_t j = 0; j < positions.rows(); j++) { 		Eb += positions.row(i) *
  //A * positions.transpose().col(j);
  //	}
  //}

  double bE;
  double sE;
  double pE;
  double cE = 0;
  double totalE;

  if (f.mesh.hasBoundary()) {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = f.H - f.H0;
    double A_difference = f.surfaceArea - f.targetSurfaceArea;
    double V_difference = f.volume - f.refVolume * f.P.Vt;

    bE = (f.P.Kb * f.M * (H_difference.array() * H_difference.array()).matrix())
             .sum();
    sE = f.P.Ksg * A_difference;
    pE = f.P.Kv * V_difference;

    if (f.isProtein) {
      double cE = (f.M * f.P.epsilon * f.proteinDensity.raw()).sum();
    }

    totalE = bE + sE + pE + cE;

  } else {

    Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = f.H - f.H0;
    double A_difference = f.surfaceArea - f.targetSurfaceArea;
    double V_difference = f.volume - f.refVolume * f.P.Vt;

    bE = (f.P.Kb * f.M * (H_difference.array() * H_difference.array()).matrix())
              .sum();
    sE = f.P.Ksg * A_difference * A_difference / f.targetSurfaceArea / 2;
    pE = f.P.Kv * V_difference * V_difference / (f.refVolume * f.P.Vt) / 2;

    if (f.isProtein) {
      double cE = (f.M * f.P.epsilon * f.proteinDensity.raw()).sum(); 
    }

    totalE = bE + sE + pE + cE;
  }

  std::tuple<double, double, double, double, double> output(totalE, bE, sE, pE, cE);
  
  return output;
}
} // namespace integration
} // namespace ddgsolver
