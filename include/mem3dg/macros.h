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

/**
 * @file macros.h
 * @brief Helpful macros for logging and setting visibility
 */

#pragma once

#include <exception>
#include <iostream>
#include <sstream>
// #include <stacktrace>  # Introduced in C++ 23!

#ifdef MEM3DG_WITH_LIBUNWIND
#define UNW_LOCAL_ONLY
#include <libunwind.h>

inline void show_backtrace() {
  unw_cursor_t cursor;
  unw_context_t context;

  // Initialize cursor to current frame for local unwinding.
  unw_getcontext(&context);
  unw_init_local(&cursor, &context);

  // Unwind frames one by one, going up the frame stack.
  while (unw_step(&cursor) > 0) {
    unw_word_t offset, ip;
    // get instruction pointer
    unw_get_reg(&cursor, UNW_REG_IP, &ip);
    if (ip == 0) {
      break;
    }
    printf("0x%lx:", ip);

    char sym[256];
    if (unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
      printf(" (%s+0x%lx)\n", sym, offset);
    } else {
      printf(" -- error: unable to obtain symbol name for this frame\n");
    }
  }
}
#else
#define show_backtrace() ((void)0)
#endif

/// Mem3dg namespace
namespace mem3dg {
/// @brief  detail namespace
namespace detail {

/**
 * @brief Construct text including file, line, and message then throw
 * runtime_error
 *
 * @tparam T pack of typenames for objects to print
 * @param function name of function
 * @param file name of file
 * @param line current line number of offense
 * @param ts pack of objects to print
 */
template <typename... T>
void throw_runtime_error(const char *function, const char *file, const int line,
                         T &&...ts) {
  std::stringstream ss;
  ss << "Error: ";
  [[maybe_unused]] int dummy[] = {0, ((ss << std::forward<T>(ts)), 0)...};
  ss << " in function " << function << " at " << file << ":" << line;
  show_backtrace();
  throw std::runtime_error(ss.str());
}

/**
 * @brief Construct text including file, line, and message then print to stderr
 *
 * @tparam T pack of typenames for objects to print
 * @param function name of function
 * @param file name of file
 * @param line current line number of offense
 * @param ts pack of objects to print
 */
template <typename... T>
void runtime_warning(const char *function, const char *file, const int line,
                     T &&...ts) {
  std::stringstream ss;
  ss << "Warning: ";
  [[maybe_unused]] int dummy[] = {0, ((ss << std::forward<T>(ts)), 0)...};
  ss << " in function " << function << " at " << file << ":" << line;
  std::cerr << ss.str() << std::endl;
}

/**
 * @brief Construct text including file, line, and message then print to stdout
 *
 * @tparam T pack of typenames for objects to print
 * @param function name of function
 * @param file name of file
 * @param line current line number of offense
 * @param ts pack of objects to print
 */
template <typename... T>
void runtime_message(const char *function, const char *file, const int line,
                     T &&...ts) {
  std::stringstream ss;
  ss << "Message: ";
  [[maybe_unused]] int dummy[] = {0, ((ss << std::forward<T>(ts)), 0)...};
  static_cast<void>(dummy); // Avoid warning for unused variable
  ss << " in function " << function << " at " << file << ":" << line;
  std::cout << ss.str() << std::endl;
}

} // namespace detail
} // namespace mem3dg

#ifdef _MSC_VER
/// @brief Alias for msvc compatibility
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

/// @brief Helper for raising runtime errors
#define mem3dg_runtime_error(...)                                              \
  mem3dg::detail::throw_runtime_error(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                                      __VA_ARGS__);

/// @brief Helper for printing runtime warnings
#define mem3dg_runtime_warning(...)                                            \
  mem3dg::detail::runtime_warning(__PRETTY_FUNCTION__, __FILE__, __LINE__,     \
                                  __VA_ARGS__);

#ifdef NO_MEM3DG_SAFETY_CHECKS
#define MEM3DG_SAFETY_ASSERT(CONDITION, MSG)
#else
#define MEM3DG_SAFETY_ASSERT(CONDITION, MSG)                                   \
  if (!(CONDITION)) {                                                          \
    throw std::runtime_error("MEM3DG_SAFETY_ASSERT FAILURE from " +            \
                             std::string(__FILE__) + ":" +                     \
                             std::to_string(__LINE__) + " - " + (MSG));        \
  }
#endif

#ifdef NDEBUG
#define mem3dg_debug_message(...) ((void)0)
#else
/**
 * @brief Helper for printing debug messages
 *
 * If NDEBUG is defined then this function is void.
 */
#define mem3dg_debug_message(...)                                              \
  mem3dg::detail::runtime_message(__PRETTY_FUNCTION__, __FILE__, __LINE__,     \
                                  __VA_ARGS__);
#endif

namespace detail {
/**
 * @brief Terminating case of helpful print
 *
 * @tparam T typename of argument to print
 * @param ss reference to active ostringstream
 * @param t object to print
 */
template <typename T> void mem3dg_print(std::ostringstream &ss, T &&t) {
  ss << std::forward<T>(t) << " ";
}

/**
 * @brief Recursive condition of helpful print
 *
 * @tparam T typename of current argument
 * @tparam TS pack of subsequent types
 * @param ss reference to active ostringstream
 * @param t current object to output
 * @param ts pack of following objects
 */
template <typename T, typename... TS>
void mem3dg_print(std::ostringstream &ss, T &&t, TS &&...ts) {
  ss << std::forward<T>(t) << " ";
  mem3dg_print(ss, std::forward<TS>(ts)...);
}
} // namespace detail

/**
 * @brief A helpful print function which will auto add spaces between printed
 * objects.
 *
 * @tparam TS pack of types for objects to print
 * @param ts pack of objects
 */
template <typename... TS> void mem3dg_print(TS &&...ts) {
  std::ostringstream ss;
  detail::mem3dg_print(ss, std::forward<TS>(ts)...);
  std::cout << ss.str() << std::endl;
}

/**
 * @def DLL_PUBLIC
 * Sets public visibility for function
 */
/**
 * @def DLL_LOCAL
 * Sets hidden visibility for function
 */

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
