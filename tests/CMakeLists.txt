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

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG origin/main
  GIT_SHALLOW TRUE
)
# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
endif()
add_subdirectory(
  ${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} SYSTEM EXCLUDE_FROM_ALL
)

# FetchContent_MakeAvailable(googletest)

# Build the tests
set(MEM3DG_TEST_SRCS src/main_test.cpp src/product_test.cpp src/force_test.cpp
                     src/integrator_test.cpp src/mutable_trajfile_test.cpp
)

add_executable(Mem3DG-tests "${MEM3DG_TEST_SRCS}")
target_link_libraries(Mem3DG-tests mem3dg gtest_main)

add_test(NAME Mem3DG_Main_Tests COMMAND Mem3DG-tests)

# Configure testing of Python module
find_package(pytest)
if(NOT PYTEST_FOUND AND BUILD_PYMEM3DG)
  message(
    WARNING
      "pytest could not be found! Please install pytest to run python tests."
  )
elseif(BUILD_PYMEM3DG)
  # Add a test for pytest
  add_test(
    NAME Mem3DG_Python_Tests
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY $<TARGET_FILE_DIR:_core>
  )
endif()

# ##############################################################################
# Main DDG target
# ##############################################################################
add_executable(ddg scratch.cpp)
target_link_libraries(ddg PRIVATE mem3dg)

add_executable(ddg_patch scratch_patch.cpp)
target_link_libraries(ddg_patch PRIVATE mem3dg)

add_executable(scratch_ctl scratch_ctl.cpp)
target_link_libraries(scratch_ctl PRIVATE mem3dg)
