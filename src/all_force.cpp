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

#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "mem3dg/solver/system.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::getAllForces() {
  if (P.Kb != 0) {
    getBendingPressure();
  }
  if (P.Kv != 0) {
    getInsidePressure();
  }
  if (P.Ksg != 0) {
    getCapillaryPressure();
  }
  if (P.eta != 0) {
    getLineTensionPressure();
  }
  if ((P.Kse != 0) || (P.Ksl != 0) || (P.Kst != 0)) {
    getRegularizationForce();
  }
  if ((P.gamma != 0) || (P.sigma != 0)){
    getDPDForces();
  }
  if (isProtein){
    getChemicalPotential();
  }
  if (P.Kf != 0){
    getExternalPressure();
  }
} // namespace ddgsolver
} // end namespace ddgsolver
