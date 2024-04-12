#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2020- The Mem3DG Authors
# and the project initiators Cuncheng Zhu, Christopher T. Lee, and
# Padmini Rangamani.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Please help us support Mem3DG development by citing the research
# papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
# for more information.

# Use CMAKE_MSVC_RUNTIME_LIBRARY to control linked runtime
cmake_policy(SET CMP0091 NEW)
# Honor visibility for all target types
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW) # Force globally
cmake_policy(SET CMP0063 NEW)
# respect project version
cmake_policy(SET CMP0048 NEW)

# 3.25 required for SYSTEM target property
cmake_minimum_required(VERSION 3.25)

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(DefaultBuildTypeAndFlags)
include(Macros)

get_version_from_git()
project(
  Mem3DG
  VERSION ${VERSION_RELEASE}
  DESCRIPTION "Membrane Dynamics in 3D using Discrete Differential Geometry"
  HOMEPAGE_URL "https://rangamanilabucsd.github.io/Mem3DG/"
  LANGUAGES CXX
)
message(STATUS "Mem3DG version: ${VERSION}")

# Configure version.cpp to give access to version in code
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
set(VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Require c++20 and standard libraries
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize output directory locations
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Check if DDG is being used directly or via add_subdirectory
set(DDG_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(DDG_MASTER_PROJECT ON)
endif()

print_debug_messages(CXX)

include(FetchContent)

option(BUILD_PYMEM3DG "Build the python extensions?" ON)
option(WITH_NETCDF "Build with NetCDF (binary trajectory output)?" ON)
option(BUILD_MEM3DG_DOCS "Configure documentation" OFF)
option(M3DG_GET_OWN_EIGEN "Download own Eigen" ON)
option(M3DG_GET_OWN_PYBIND11 "Download own pybind11" ON)
option(WITH_LIBUNWIND "Link libunwind for stack traces" OFF)
option(WITH_OPENMP "Build with OpenMP support" OFF)
option(LINK_PROFILER "Link profiler gperftools" OFF)

option(MEM3DG_PEDANTIC "Be extremely pedantic while compiling" OFF)

mark_as_advanced(FORCE MEM3DG_PEDANTIC)

# ##############################################################################
# BUNDLED LIBRARIES & EXTERNAL LIBS
# ##############################################################################
add_subdirectory(libraries)

set(LINKED_LIBS geometry-central pcg::pcg igl::core)

if(WITH_NETCDF)
  if(VCPKG_TARGET_TRIPLET)
    message(STATUS "VCPKG Triplet: ${VCPKG_TARGET_TRIPLET}")
  endif()
  find_package(netCDF "4.3.1" REQUIRED)
  message(DEBUG "netCDF version: ${netCDF_VERSION}")
  message(DEBUG "netCDF dir: ${netCDF_DIR}")
  message(DEBUG "netCDF include dir: ${netCDF_INCLUDE_DIR}")
  message(DEBUG "netCDF libraries: ${netCDF_LIBRARIES}")
  list(APPEND LINKED_LIBS NetCDF::NetCDF)

  find_package(netcdf-cxx4 REQUIRED)
  message(DEBUG "netcdf-cxx4 include: ${netcdf-cxx4_INCLUDE_DIRS}")
  message(DEBUG "netcdf-cxx4 library: ${netcdf-cxx4_LIBRARIES}")
  message(DEBUG "netcdf-cxx4 version: ${netcdf-cxx4_VERSION}")
  list(APPEND LINKED_LIBS NetCDF::NetCDF-cxx4)
endif()

# ##############################################################################
# DDG SOLVER LIBRARY
# ##############################################################################
add_subdirectory(include) # Get list of headers
add_subdirectory(src) # Get list of sources
list(APPEND DDG_SOURCES ${VERSION_FILE})

# mem3dg object library
add_library(mem3dg_objlib OBJECT ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(
  mem3dg_objlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  # list(APPEND LINKED_LIBS OpenMP::OpenMP_CXX)
  target_link_libraries(
    mem3dg_objlib PRIVATE $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>
  )
endif()

if(LINK_PROFILER)
  find_package(Gperftools REQUIRED)
  list(APPEND LINKED_LIBS ${GPERFTOOLS_PROFILER})
  target_include_directories(mem3dg_objlib PUBLIC ${GPERFTOOLS_INCLUDE_DIR})
  target_compile_definitions(mem3dg_objlib PUBLIC -DMEM3DG_WITH_GPERFTOOLS)
endif()

target_link_libraries(mem3dg_objlib PUBLIC ${LINKED_LIBS})
set_target_properties(mem3dg_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(WITH_LIBUNWIND)
  target_compile_definitions(mem3dg_objlib PUBLIC -DMEM3DG_WITH_LIBUNWIND)
endif()

if(WITH_NETCDF)
  target_compile_definitions(mem3dg_objlib PUBLIC -DMEM3DG_WITH_NETCDF)
endif()

if(MEM3DG_PEDANTIC)
  if(NOT MSVC)
    target_compile_options(
      mem3dg_objlib PRIVATE "-Wall" "-Wextra" "-Wshadow" "-Wnon-virtual-dtor"
                            "-pedantic"
    )
  endif(NOT MSVC)
endif()

# mem3dg library
add_library(mem3dg SHARED $<TARGET_OBJECTS:mem3dg_objlib>)
target_link_libraries(mem3dg PUBLIC mem3dg_objlib)

if(BUILD_PYMEM3DG)
  find_package(
    Python
    COMPONENTS Interpreter Development
    REQUIRED
  )
  add_subdirectory(python_src)
endif()

enable_testing()
add_subdirectory(tests)

if(BUILD_MEM3DG_DOCS)
  add_subdirectory(docs)
endif()

include(GNUInstallDirs)
# Install targets and headers
install(
  TARGETS mem3dg
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
  DIRECTORY include/mem3dg
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.h"
)
