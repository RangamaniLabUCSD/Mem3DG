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

set(DDG_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/typetraits.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/macros.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/type_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/meshops.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/version.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/constants.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/mesh_io.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/system.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/geometry.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/parameters.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/forces.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/mesh_process.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/trajfile_constants.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/trajfile.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/mutable_trajfile.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/integrator/integrator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/integrator/forward_euler.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/integrator/conjugate_gradient.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mem3dg/solver/integrator/velocity_verlet.h"
    PARENT_SCOPE
)
