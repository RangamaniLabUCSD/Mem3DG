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

function(eigen3checker EIGEN_LOCATION EIGEN3_FIND_VERSION)
  # Function looks for Eigen3 at EIGEN_LOCATION Eigen version much be >=
  # EIGEN3_FIND_VERSION
  #
  # If Eigen is found it sets the following variables on PARENT_SCOPE
  #
  # EIGEN3_FOUND       - system has eigen lib with correct version
  # EIGEN3_INCLUDE_DIR - the eigen include directory EIGEN3_VERSION     - eigen
  # version
  #

  set(EIGEN3_FOUND
      false
      PARENT_SCOPE
  )

  # Search for the signature_of_eigen3_matrix_library

  # Search first for just the requested path
  find_path(
    EIGEN3_INCLUDE_DIR
    NAMES signature_of_eigen3_matrix_library
    PATHS ${EIGEN_LOCATION}
    NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
  )

  # Now search more broadly (will do nothing if the search above succeeded)
  find_path(
    EIGEN3_INCLUDE_DIR
    NAMES signature_of_eigen3_matrix_library
    PATH_SUFFIXES eigen3 eigen
  )

  if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")

    # Parse version from Macros.h
    file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h"
         _eigen3_version_header
    )

    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)"
                 _eigen3_world_version_match "${_eigen3_version_header}"
    )
    set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)"
                 _eigen3_major_version_match "${_eigen3_version_header}"
    )
    set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)"
                 _eigen3_minor_version_match "${_eigen3_version_header}"
    )
    set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")

    set(EIGEN3_VERSION
        ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION}
    )

    # message(STATUS "Eigen Version: ${EIGEN3_VERSION}")
    if(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
      # Complain about the inadequate version
      message(
        STATUS
          "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, "
          "but at least version ${EIGEN3_FIND_VERSION} is required"
      )
    else(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
      # Suitable version of Eigen3 found. Set variables
      set(EIGEN3_VERSION
          ${EIGEN3_VERSION}
          PARENT_SCOPE
      )
      set(EIGEN3_FOUND
          true
          PARENT_SCOPE
      )
      set(EIGEN3_INCLUDE_DIR
          ${EIGEN3_INCLUDE_DIR}
          PARENT_SCOPE
      )
    endif(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
  else()
    # Purge EIGEN3_INCLUDE_DIR value since it's invalid Required so find_path
    # will actually search instead of defaulting
    set(EIGEN3_INCLUDE_DIR
        "EIGEN3_INCLUDE_DIR-NOTFOUND"
        PARENT_SCOPE
    )
  endif()
endfunction()
