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

#include "ddgsolver/force.h"
#include "ddgsolver/integrator.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>

namespace ddgsolver {
namespace integration {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

double getBendingEnergy(Force &f) {
  // comment: this may not be useful, the convergence can be be tested by
  // checking its derivative which is the forces excluding the DPD forces. The
  // energy trajectory of could actually numerically integrated by post
  // processing after saving all forces during the iterations.

  // auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);
  // auto A = f.L.transpose() * f.M_inv * f.L; //could further cached since the
  // same for all time pt for (size_t i = 0; i < positions.rows(); i++) {
  // for
  //(size_t j = 0; j < positions.rows(); j++) { 		Eb += positions.row(i) *
  //A * positions.transpose().col(j);
  //	}
  //}

  auto difference = f.M_inv * f.H - f.H0;
  double bE =
      (f.P.Kb * f.M * (difference.array() * difference.array()).matrix()).sum();
  /// stretching energy
  // double sE =

  return bE;
}
} // namespace integration
} // namespace ddgsolver
