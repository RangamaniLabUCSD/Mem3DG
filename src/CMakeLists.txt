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

set(DDG_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/icosphere.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/bending_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/stretching_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pressure_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dpd_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/stormer_verlet.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/velocity_verlet.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/bending_energy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/external_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/vesicle_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/patch_force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/log_file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/trajfile.cpp"
    PARENT_SCOPE
)
