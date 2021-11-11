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

#pragma once

// #include <cassert>
// #include <geometrycentral/surface/surface_mesh.h>
// #include <geometrycentral/surface/vertex_position_geometry.h>
// #include <geometrycentral/utilities/eigen_interop_helpers.h>

// #include <Eigen/Core>
// #include <Eigen/SparseLU>

// #include <functional>
// #include <math.h>
// #include <vector>
// #include <iostream>

// #include "geometrycentral/surface/halfedge_element_types.h"
// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/utilities/vector2.h"
// #include "geometrycentral/utilities/vector3.h"

#include "mem3dg/solver/physics/physic.h"

// namespace gc = ::geometrycentral;
// namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {

class DLL_PUBLIC StretchingPhysic : public Physic {
public:
  StretchingPhysic(){};

  inline void compute(VariablePhysicsSystem& S) override {
    std::cout << "Computing Stretching" << std::endl;
  };
  
  inline void initialize(VariablePhysicsSystem& S) override {
    std::cout << "Initialize Stretching" << std::endl;
  };
};
} // namespace solver
} // namespace mem3dg
