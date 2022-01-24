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

#include <exception>
#include <sstream>

namespace mem3dg {
namespace internal {
template <typename... T>
void throw_runtime_error(const char *function, const char *file, const int line,
                         T &&...ts) {
  std::stringstream ss;
  ss << "Error: ";
  int dummy[] = {0, ((ss << std::forward<T>(ts)), 0)...};
  static_cast<void>(dummy); // Avoid warning for unused variable
  ss << " in function " << function << " at " << file << ":" << line;
  throw std::runtime_error(ss.str());
}
template <typename... T>
void throw_runtime_message(const char *function, const char *file, const int line,
                         T &&...ts) {
  std::stringstream ss;
  ss << "Message: ";
  int dummy[] = {0, ((ss << std::forward<T>(ts)), 0)...};
  static_cast<void>(dummy); // Avoid warning for unused variable
  ss << " in function " << function << " at " << file << ":" << line;
  std::cout << ss.str() << std::endl;
}
} // namespace internal
} // namespace mem3dg

#ifdef _MSC_VER 
  #define __PRETTY_FUNCTION__ __FUNCSIG__ 
#endif

#define mem3dg_runtime_error(...)                                              \
  mem3dg::internal::throw_runtime_error(__PRETTY_FUNCTION__, __FILE__,         \
                                        __LINE__, __VA_ARGS__);

#define mem3dg_runtime_message(...)                                              \
  mem3dg::internal::throw_runtime_message(__PRETTY_FUNCTION__, __FILE__,         \
                                        __LINE__, __VA_ARGS__);

#if defined _WIN32 || defined __CYGWIN__
#ifdef _DLL
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__((dllexport))
#else
#define DLL_PUBLIC                                                             \
  __declspec(                                                                  \
      dllexport) // Note: actually gcc seems to also supports this syntax.
#endif
#else
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__((dllimport))
#else
#define DLL_PUBLIC                                                             \
  __declspec(                                                                  \
      dllimport) // Note: actually gcc seems to also supports this syntax.
#endif
#endif
#define DLL_LOCAL
#else
#if __GNUC__ >= 4
#define DLL_PUBLIC __attribute__((visibility("default")))
#define DLL_LOCAL __attribute__((visibility("hidden")))
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif
#endif
