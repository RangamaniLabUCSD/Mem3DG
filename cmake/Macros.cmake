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

#[[
Macro sets the following variables
  * `VERSION`
  * `VERSION_RELEASE`
  * `VERSION_PRE`
  * `VERSION_DEV`
  * `VERSION_SHA1`
  * `VERSION_DIRTY`
#]]
macro(get_version_from_git)
  # Look for version from GIT
  include(GetGitRevisionDescription)
  git_describe_working_tree(VERSION --tags --always --long)

  if(VERSION MATCHES [[.*NOTFOUND]])
    # else() Check for PKG-INFO from sdist
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/PKG-INFO)
      file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/PKG-INFO VERSION
           REGEX "^Version: (.*)$"
      )
      string(REPLACE "Version: " "" VERSION ${VERSION})
      set(USING_PKG_INFO 1)
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git_archival.txt)
      file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/.git_archival.txt VERSION
           REGEX "^describe-name: (.*)$"
      )
      string(REPLACE "describe-name: " "" VERSION ${VERSION})
    endif()
  endif()
  # message(STATUS ${VERSION})
  string(
    REGEX
      MATCH
      [[^v?([0-9]+(\.[0-9]+)*)[_.-]?(rc[0-9]*|a[0-9]*|b[0-9]*|c[0-9]*)?[_.-]?(dev)?([0-9]+)?[_.+-]?([a-z0-9]+)?[_.-]?(d[0-9]+)?[-+]?(dirty)?$]]
      MATCH_RESULT
      "${VERSION}"
  )
  # cmake-format: off
  # -- MATCH 0: 0.0.6rc3.post1.dev13+gc0963e4.d20240124
  # -- MATCH 1: 0.0.6  # Release version
  # -- MATCH 2: .6     # Last of release
  # -- MATCH 3: rc3    # Pre-release
  # -- MATCH 4: dev   # dummy
  # -- MATCH 5: 13    # Dev version
  ## Local version
  # -- MATCH 6: gc0963e4  # SHA1
  # -- MATCH 7: d20240124
  # -- MATCH 8: dirty

  if(MATCH_RESULT)
    # foreach(_TMP RANGE 9)
    #   message(STATUS "MATCH ${_TMP}: ${CMAKE_MATCH_${_TMP}}")
    # endforeach()
    # cmake-format: on

    set(VERSION_RELEASE ${CMAKE_MATCH_1})
    set(VERSION_PRE ${CMAKE_MATCH_3})
    set(VERSION_DEV ${CMAKE_MATCH_5})
    set(VERSION_SHA1 ${CMAKE_MATCH_6})
    set(VERSION_DIRTY ${CMAKE_MATCH_8})

    string(REGEX MATCHALL [[\.?[0-9]+]] VERSION_RELEASE_L ${VERSION_RELEASE})
    foreach(X ${VERSION_RELEASE_L})
      string(REPLACE "." "" X ${X})
      list(APPEND VERSION_RELEASE_L_tmp ${X})
    endforeach()
    set(VERSION_RELEASE_L ${VERSION_RELEASE_L_tmp})

    # Split out prerelease version
    if(VERSION_PRE)
      string(REGEX MATCH "([a-z]*)([0-9]+)" MATCH_RESULT "${VERSION_PRE}")
      set(VERSION_PRE_L ${CMAKE_MATCH_1})
      set(VERSION_PRE_N ${CMAKE_MATCH_2})
    endif()
    # message(STATUS "${VERSION_LEAST} ${VERSION_PRE_L} ${VERSION_PRE_N}")
  else()
    message(WARNING "Version could not be determined")
    set(VERSION_RELEASE "0.0.0")
    set(VERSION_PRE "")
    set(VERSION_DEV "")
    set(VERSION_SHA1 "")
    set(VERSION_DIRTY "")
  endif()

  # If dirty increment version
  if((VERSION_DIRTY OR VERSION_DEV) AND NOT USING_PKG_INFO)
    if(VERSION_PRE)
      math(EXPR VERSION_PRE_N "${VERSION_PRE_N} + 1")
      set(VERSION_PRE "${VERSION_PRE_L}${VERSION_PRE_N}")
    else()
      list(POP_BACK VERSION_RELEASE_L X)
      math(EXPR X "${X} + 1")
      list(APPEND VERSION_RELEASE_L ${X})
      string(JOIN "." VERSION_RELEASE ${VERSION_RELEASE_L})
    endif()
  endif()

  # Reconstruct version
  set(VERSION ${VERSION_RELEASE})
  set(VERSION_SHORT ${VERSION_RELEASE})
  if(VERSION_PRE)
    string(APPEND VERSION ${VERSION_PRE})
    string(APPEND VERSION_SHORT ${VERSION_PRE})
  endif()
  if(VERSION_DIRTY
     OR VERSION_DEV
     OR USING_PKG_INFO
  )
    if(NOT VERSION_DEV STREQUAL "")
      string(APPEND VERSION ".dev${VERSION_DEV}")
    endif()
    if(VERSION_SHA1)
      string(APPEND VERSION "+${VERSION_SHA1}")
    endif()
  endif()
endmacro(get_version_from_git)

macro(print_debug_messages LANG)
  message(DEBUG "CMAKE_${LANG}_FLAGS is: ${CMAKE_${LANG}_FLAGS}")
  message(DEBUG "CMAKE_${LANG}_FLAGS_DEBUG is: ${CMAKE_${LANG}_FLAGS_DEBUG}")
  message(DEBUG
          "CMAKE_${LANG}_FLAGS_RELEASE is: ${CMAKE_${LANG}_FLAGS_RELEASE}"
  )
  message(
    DEBUG
    "CMAKE_${LANG}_FLAGS_RELWITHDEBINFO is: ${CMAKE_${LANG}_FLAGS_RELWITHDEBINFO}"
  )
  message(DEBUG
          "CMAKE_${LANG}_FLAGS_MINSIZEREL is: ${CMAKE_${LANG}_FLAGS_MINSIZEREL}"
  )
  message(DEBUG "Build type: ${CMAKE_BUILD_TYPE}")
  message(DEBUG "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE})
  message(DEBUG "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS})
endmacro(print_debug_messages)
