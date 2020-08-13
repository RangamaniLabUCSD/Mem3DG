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

# ##############################################################################
# GET PYBIND11
# ##############################################################################
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.5.0
  GIT_SHALLOW TRUE
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src" BINARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
  FetchContent_Populate(pybind11)
  add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

# mem3dg static library
add_library(
  mem3dg_static STATIC EXCLUDE_FROM_ALL $<TARGET_OBJECTS:mem3dg_objlib>
)
target_link_libraries(mem3dg_static PUBLIC mem3dg_objlib)
set_target_properties(
  mem3dg_static PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET
                                                        hidden
)

set(PYMEM3DG_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/pyddg.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ddgsolver.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization.cpp"
)

pybind11_add_module(pymem3dg MODULE ${PYMEM3DG_SOURCES})
target_link_libraries(pymem3dg PRIVATE mem3dg_static polyscope)


file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
             INPUT ${CMAKE_SOURCE_DIR}/cmake/__init__.py.in)

if(SKBUILD)
  install(TARGETS pymem3dg LIBRARY DESTINATION pymem3dg)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION pymem3dg)
else()
  message(
    WARNING
      "CMake install no longer installs Pymem3DG to your Python site-packages. "
      "Please configure and install using 'python setup.py install' instead."
  )
endif()
