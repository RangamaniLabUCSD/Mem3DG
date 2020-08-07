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

#include "geometrycentral/utilities/vector3.h"
#include "macros.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;

/**
 * @brief Construct an icosphere mesh in PolygonSoup form
 *
 * @param coords    Reference to vector of vertex coordinates
 * @param polygons  Polygons of the mesh
 * @param n         Iterations of quadrisections to perform
 */
DLL_PUBLIC void icosphere(std::vector<gc::Vector3> &coords,
                          std::vector<std::vector<std::size_t>> &polygons,
                          int n, double R = 1);

DLL_PUBLIC void tetrahedron(std::vector<gc::Vector3> &coords,
                            std::vector<std::vector<std::size_t>> &polygons);
} // end namespace ddgsolver
