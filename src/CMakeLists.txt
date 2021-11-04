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
    "${CMAKE_CURRENT_SOURCE_DIR}/mesh_io.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/visualization.cpp"
    
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/init.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/force.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/energy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/parameters.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/regularization.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/mesh_process.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/trajfile.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/mutable_trajfile.cpp"

    "${CMAKE_CURRENT_SOURCE_DIR}/solver/integrator/integrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/integrator/BFGS.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/integrator/velocity_verlet.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/integrator/forward_euler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/integrator/conjugate_gradient.cpp"
    PARENT_SCOPE
)
