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

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.10.0
  GIT_SHALLOW TRUE
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-src" 
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
)
# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  # Add googletest directly to our build. This defines the gtest and gtest_main
  # targets.
  add_subdirectory(
    ${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL
  )
endif()

set_property(TARGET gtest PROPERTY CXX_STANDARD 11)

# Build the tests
set(TEST_SRCS src/main_test.cpp src/dotproduct_test.cpp
              src/force_calculation_test.cpp
)

add_executable(Mem3DG-tests "${TEST_SRCS}")
target_link_libraries(Mem3DG-tests mem3dg gtest_main)

add_test(NAME Mem3DG_Main_Tests COMMAND Mem3DG-tests)


# ##############################################################################
# Main DDG target
# ##############################################################################
# add_executable(ddg scratch.cpp) target_link_libraries(ddg PRIVATE
# mem3dg polyscope) set_target_properties(ddg PROPERTIES
# CXX_VISIBILITY_PRESET hidden)
