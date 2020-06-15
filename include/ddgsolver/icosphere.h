

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
               std::vector<std::vector<std::size_t>> &polygons, int n);
} // end namespace ddgsolver
