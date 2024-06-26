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

# ##############################################################################
# GET PYBIND11
# ##############################################################################
if(M3DG_GET_OWN_PYBIND11)
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.10.4
    GIT_SHALLOW TRUE
    SYSTEM SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src" BINARY_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
  )
  FetchContent_MakeAvailable(pybind11)
else()
  find_package(pybind11 REQUIRED CONFIG)
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
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_forces.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_integrators.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_geometry.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_mutators.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_system.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_parameters.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pymem3dg_energy.cpp"
)

pybind11_add_module(_core MODULE ${PYMEM3DG_SOURCES})
target_link_libraries(_core PRIVATE mem3dg_static)

if(WITH_NETCDF)
  target_compile_definitions(_core PUBLIC -DMEM3DG_WITH_NETCDF)
endif()

install(TARGETS _core LIBRARY DESTINATION pymem3dg/ COMPONENT PythonModule)

# install(TARGETS _core DESTINATION .) file(GENERATE OUTPUT
# ${CMAKE_CURRENT_BINARY_DIR}/__init__.py INPUT
# ${CMAKE_SOURCE_DIR}/cmake/__init__.py.in)

# if(SKBUILD) install(TARGETS pymem3dg LIBRARY DESTINATION pymem3dg)
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION pymem3dg)
# endif()
