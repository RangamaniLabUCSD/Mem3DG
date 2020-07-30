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

project(
  DDG
  VERSION 0.0.1
  LANGUAGES CXX
)
if(MSVC)
    set(CMAKE_C_FLAGS_DEBUG_INIT "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_C_FLAGS_MINSIZEREL_INIT     "/MT /O1 /Ob1 /D NDEBUG")
    set(CMAKE_C_FLAGS_RELEASE_INIT        "/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "/MT /Zi /O2 /Ob1 /D NDEBUG")
else()
    set(CMAKE_C_FLAGS_RELEASE_INIT        "-O2 -DNDEBUG")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/MT /O1 /Ob1 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/MT /Zi /O2 /Ob1 /D NDEBUG")
else()
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "-O2 -DNDEBUG")
endif()

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(DefaultBuildType)

# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize output directory locations
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Required to tell MSVC to export symbols
# if(MSVC)
#  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
# endif()

# Check if DDG is being used directly or via add_subdirectory
set(DDG_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(DDG_MASTER_PROJECT ON)
endif()

# find_package(Git REQUIRED)
include(FetchContent)

option(BUILD_PYDDG "Build the python extensions?" ON)

# set(CMAKE_CXX_VISIBILITY_PRESET hidden)

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
# BUNDLED LIBRARIES
# ##############################################################################
add_subdirectory(libraries)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/conservative_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/log_file.cpp

)

# ddgsolver object library
add_library(ddgsolver_objlib OBJECT ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(
  ddgsolver_objlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ddgsolver_objlib PUBLIC geometry-central pcg::pcg)
set_target_properties(ddgsolver_objlib PROPERTIES POSIION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden) # VISIBILITY_INLINES_HIDDEN ON)

# ddgsolver library
add_library(ddgsolver SHARED $<TARGET_OBJECTS:ddgsolver_objlib>)
target_link_libraries(ddgsolver PUBLIC ddgsolver_objlib)

if(BUILD_PYDDG)
  add_subdirectory(pyddg)
endif()

# ##############################################################################
# Main DDG target
# ##############################################################################
enable_testing()
add_subdirectory(tests)
