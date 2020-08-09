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

set(DDG_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/constants.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/ddgsolver.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/force.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/integrator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/icosphere.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/trajfile.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/typetraits.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/util.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/macros.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/meshops.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/version.h"
    PARENT_SCOPE)
