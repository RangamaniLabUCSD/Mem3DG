# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
#

cmake_minimum_required(VERSION 3.11)

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(DefaultBuildTypeAndFlags)
include(Macros)

get_version_from_git()
project(
  Mem3DG
  VERSION ${VERSION_SHORT}
  LANGUAGES CXX
)
message(STATUS "Mem3DG version: ${VERSION}")

if(SKBUILD)
    message(STATUS "The project is being built using scikit-build")
endif()

# Configure version.cpp to give access to version in code
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
set(VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
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

print_debug_messages()

include(FetchContent)

option(BUILD_PYMEM3DG "Build the python extensions?" ON)
option(WITH_NETCDF "Build with NetCDF (binary trajectory output)?" ON)
option(BUILD_MEM3DG_DOCS "Build the python extensions?" ON)

# ##############################################################################
# BUNDLED LIBRARIES & EXTERNAL LIBS
# ##############################################################################
add_subdirectory(libraries)

set(LINKED_LIBS geometry-central pcg::pcg)

if(WITH_NETCDF)
  find_package(netCDF REQUIRED)
  message(STATUS ${netCDF_VERSION})
  message(STATUS ${netCDF_DIR})
  message(STATUS ${netCDF_LIBRARIES})
  list(APPEND LINKED_LIBS NetCDF::NetCDF)

  find_package(NetCDF_cxx4)

  if(NetCDF_cxx4_FOUND)
    list(APPEND LINKED_LIBS NetCDF::NetCDF_cxx4)
  else()
    find_path(netcdf-cxx4_INCLUDE_DIRS NAMES netcdf REQUIRED)
    include(SelectLibraryConfigurations)
    find_library(netcdf-cxx4_LIBRARY_RELEASE NAMES netcdf-cxx4)
    find_library(netcdf-cxx4_LIBRARY_DEBUG NAMES netcdf-cxx4 PATH_SUFFIXES debug/lib)
    select_library_configurations(netcdf-cxx4)
    
    message(STATUS ${netcdf-cxx4_INCLUDE_DIRS})
    message(STATUS ${netcdf-cxx4_LIBRARY})

    list(APPEND LINKED_LIBS ${netcdf-cxx4_LIBRARY})
  endif()
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
target_link_libraries(mem3dg_objlib PUBLIC ${LINKED_LIBS})
set_target_properties(
  mem3dg_objlib PROPERTIES POSIION_INDEPENDENT_CODE ON)
if(WITH_NETCDF)
  target_compile_definitions(mem3dg_objlib PUBLIC -DMEM3DG_WITH_NETCDF)
  target_include_directories(mem3dg_objlib PRIVATE ${netcdf-cxx4_INCLUDE_DIRS})
endif()

# mem3dg library
add_library(mem3dg SHARED $<TARGET_OBJECTS:mem3dg_objlib>)
target_link_libraries(mem3dg PUBLIC mem3dg_objlib)

if(BUILD_PYMEM3DG)
  add_subdirectory(pymem3dg)
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
