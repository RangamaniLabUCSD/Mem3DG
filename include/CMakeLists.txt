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
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/constants.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/ddgsolver.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/force.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/integrator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/icosphere.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/trajfile.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/typetraits.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/util.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/macros.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/meshops.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/ddgsolver/version.h"
    PARENT_SCOPE)
