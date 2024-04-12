#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2024- The Mem3DG Authors
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

# Tries to find Gperftools.
#
# Usage of this module as follows:
#
# find_package(Gperftools)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
# Gperftools_ROOT_DIR  Set this variable to the root installation of Gperftools
# if the module has problems finding the proper installation path.
#
# Variables defined by this module:
#
# GPERFTOOLS_FOUND              System has Gperftools libs/headers
# GPERFTOOLS_LIBRARIES          The Gperftools libraries (tcmalloc & profiler)
# GPERFTOOLS_INCLUDE_DIR        The location of Gperftools headers

find_library(
  GPERFTOOLS_TCMALLOC
  NAMES tcmalloc
  HINTS ${Gperftools_ROOT_DIR}/lib
)

find_library(
  GPERFTOOLS_PROFILER
  NAMES profiler
  HINTS ${Gperftools_ROOT_DIR}/lib
)

find_library(
  GPERFTOOLS_TCMALLOC_AND_PROFILER
  NAMES tcmalloc_and_profiler
  HINTS ${Gperftools_ROOT_DIR}/lib
)

find_path(
  GPERFTOOLS_INCLUDE_DIR
  NAMES gperftools/heap-profiler.h
  HINTS ${Gperftools_ROOT_DIR}/include
)

set(GPERFTOOLS_LIBRARIES ${GPERFTOOLS_TCMALLOC_AND_PROFILER})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Gperftools DEFAULT_MSG GPERFTOOLS_LIBRARIES GPERFTOOLS_INCLUDE_DIR
)

mark_as_advanced(
  Gperftools_ROOT_DIR GPERFTOOLS_TCMALLOC GPERFTOOLS_PROFILER
  GPERFTOOLS_TCMALLOC_AND_PROFILER GPERFTOOLS_LIBRARIES GPERFTOOLS_INCLUDE_DIR
)
