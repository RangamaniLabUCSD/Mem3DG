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
#     Padmini Rangmani (prangamani@eng.ucsd.edu)
# 

cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0063 NEW)

project(
  DDG
  VERSION 0.0.1
  LANGUAGES CXX
)

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(DefaultBuildTypeAndFlags)

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

include(FetchContent)

option(BUILD_PYDDG "Build the python extensions?" ON)
option(WITH_NETCDF "Build with NetCDF (binary trajectory output)?" OFF)

message(DEBUG "CMAKE_C_FLAGS is: ${CMAKE_C_FLAGS}")
message(DEBUG "CMAKE_C_FLAGS_DEBUG is: ${CMAKE_C_FLAGS_DEBUG}")
message(DEBUG "CMAKE_C_FLAGS_RELEASE is: ${CMAKE_C_FLAGS_RELEASE}")
message(DEBUG "CMAKE_C_FLAGS_RELWITHDEBINFO is: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
message(DEBUG "CMAKE_C_FLAGS_MINSIZEREL is: ${CMAKE_C_FLAGS_MINSIZEREL}")
message(DEBUG "CMAKE_CXX_FLAGS is: ${CMAKE_CXX_FLAGS}")
message(DEBUG "CMAKE_CXX_FLAGS_DEBUG is: ${CMAKE_CXX_FLAGS_DEBUG}")
message(DEBUG "CMAKE_CXX_FLAGS_RELEASE is: ${CMAKE_CXX_FLAGS_RELEASE}")
message(DEBUG "CMAKE_CXX_FLAGS_RELWITHDEBINFO is: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(DEBUG "CMAKE_CXX_FLAGS_MINSIZEREL is: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message(DEBUG "Build type: ${CMAKE_BUILD_TYPE}")
message(DEBUG "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE})
message(DEBUG "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS})

# ##############################################################################
# BUNDLED LIBRARIES & EXTERNAL LIBS
# ##############################################################################
add_subdirectory(libraries)

set(LINKED_LIBS geometry-central pcg::pcg)

if(WITH_NETCDF)
  find_package(NetCDF REQUIRED)
  list(APPEND LINKED_LIBS NetCDF::NetCDF)
endif()


# ##############################################################################
# DDG SOLVER LIBRARY
# ##############################################################################
add_subdirectory(include) # Get file sources

set(DDG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/icosphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bending_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stretching_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pressure_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dpd_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stormer_verlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/velocity_verlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bending_energy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/external_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vesicle_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tube_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/log_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/trajfile.cpp
)

# ddgsolver object library
add_library(ddgsolver_objlib OBJECT ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(
  ddgsolver_objlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ddgsolver_objlib PUBLIC ${LINKED_LIBS})
set_target_properties(ddgsolver_objlib PROPERTIES POSIION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden) # VISIBILITY_INLINES_HIDDEN ON)
if(WITH_NETCDF)
  target_compile_definitions(ddgsolver_objlib PUBLIC -DMEM3DG_WITH_NETCDF)
endif()

# ddgsolver library
add_library(ddgsolver SHARED $<TARGET_OBJECTS:ddgsolver_objlib>)
target_link_libraries(ddgsolver PUBLIC ddgsolver_objlib)

if(BUILD_PYDDG)
  add_subdirectory(pyddg)
endif()

# ##############################################################################
# Main DDG target
# ##############################################################################
# add_executable(ddg tests/scratch.cpp)
# target_link_libraries(ddg PRIVATE ddgsolver polyscope)
# set_target_properties(ddg PROPERTIES CXX_VISIBILITY_PRESET hidden)

enable_testing()
add_subdirectory(tests)
