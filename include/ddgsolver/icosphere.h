

#pragma once

#include "geometrycentral/utilities/vector3.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;

void icosphere(std::vector<gc::Vector3> &coords,
               std::vector<std::vector<std::size_t>> &polygons, int n);
} // end namespace ddgsolver
