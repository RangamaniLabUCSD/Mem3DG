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

#[[
Macro sets the following variables
  * `VERSION_SHORT`
  * `VERSION_INFO`
  * `VERSION_SHA1`
  * `VERSION_DIRTY`
#]]
macro(get_version_from_git)
  # Look for version from GIT
  include(GetGitRevisionDescription)
  git_describe_working_tree(VERSION)

  if(VERSION MATCHES [[.*NOTFOUND]])
    # else() Check for PKG-INFO from sdist
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/PKG-INFO)
      file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/PKG-INFO VERSION
           REGEX "^Version: (.*)$"
      )
      string(REPLACE "Version: " "" VERSION ${VERSION})
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git_archival.txt)
      file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/.git_archival.txt VERSION
           REGEX "^describe-name: (.*)$"
      )
      string(REPLACE "describe-name: " "" VERSION ${VERSION})
    endif()
  endif()
  message(STATUS ${VERSION})
  string(
    REGEX
      MATCH
      [[^v?([0-9]+(\.[0-9]+)*)[_.-]?(alpha[0-9]*|beta[0-9]*|dev[0-9]*|rc[0-9]*|a[0-9]*|b[0-9]*|c[0-9]*)?[_.-]?d?e?v?([0-9]+)?[_.+-]?([a-z0-9]+)?[_.-]?(d[0-9]+)?-?(dirty)?$]]
      MATCH_RESULT
      "${VERSION}"
  )
  if(MATCH_RESULT)
    foreach(_TMP RANGE 10)
      message(TRACE "MATCH ${_TMP}: ${CMAKE_MATCH_${_TMP}}")
    endforeach()
    set(VERSION_SHORT ${CMAKE_MATCH_1})
    set(VERSION_INFO ${CMAKE_MATCH_3})
    set(VERSION_SHA1 ${CMAKE_MATCH_5})

  else()
    message(WARNING "Version could not be determined")
    set(VERSION_SHORT "0.0.0")
    set(VERSION_INFO "")
    set(VERSION_SHA1 "")
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
