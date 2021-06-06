// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2021:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include <exception>
#include <sstream>

namespace mem3dg {

template <typename... T>
inline void throw_runtime_error(const char *file, const int line,
                          T&&... ts) {
  std::stringstream ss;
  ss << file << " on line " << line << ": ";
  int dummy[] = { 0, ((ss << std::forward<T>(ts)), 0)... };
  static_cast<void>(dummy); // Avoid warning for unused variable
  throw std::runtime_error(ss.str());
}
} // namespace mem3dg
