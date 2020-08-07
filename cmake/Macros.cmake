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

  * `VERSION_MAJOR`
  * `VERSION_MINOR`
  * `VERSION_PATCH`
  * `VERSION_INFO`
  * `VERSION_SHA1`
  * `VERSION_DIRTY`
  * `VERSION_SHORT`
#]]
macro(get_version_from_git)
  # Look for version from GIT
  include(GetGitRevisionDescription)
  git_describe(VERSION --tags --dirty=.dirty --always)

  # parse the version information into pieces.
  string(
    REGEX
      MATCH
      "^v([0-9]+)\\.([0-9]+)\\.([0-9]+)-*(alpha|beta|dev|)-*([A-Za-z0-9_-]*)\\.*(dirty|)$"
      MATCH_RESULT
      "${VERSION}"
  )

  if(MATCH_RESULT)
    # message(STATUS "Results: ${CMAKE_MATCH_0} ${CMAKE_MATCH_1}
    # ${CMAKE_MATCH_2} ${CMAKE_MATCH_3} ${CMAKE_MATCH_4} ${CMAKE_MATCH_5}
    # ${CMAKE_MATCH_6}")
    set(VERSION_MAJOR ${CMAKE_MATCH_1})
    set(VERSION_MINOR ${CMAKE_MATCH_2})
    set(VERSION_PATCH ${CMAKE_MATCH_3})
    set(VERSION_INFO ${CMAKE_MATCH_4})
    set(VERSION_SHA1 ${CMAKE_MATCH_5})
    set(VERSION_DIRTY ${CMAKE_MATCH_6})

    set(VERSION_DUMP "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
    if(VERSION_INFO)
      string(CONCAT VERSION_DUMP ${VERSION_DUMP} "-${VERSION_INFO}")
    endif()
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VERSION.in
      ${CMAKE_CURRENT_SOURCE_DIR}/VERSION
    )
  else()
    message(STATUS "No GIT VCS found pulling version from file")
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/VERSION VERSION)
    string(STRIP "${VERSION}" VERSION)
    string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)-*(alpha|beta|dev|)$"
                 MATCH_RESULT "${VERSION}"
    )
    if(MATCH_RESULT)
      set(VERSION_MAJOR ${CMAKE_MATCH_1})
      set(VERSION_MINOR ${CMAKE_MATCH_2})
      set(VERSION_PATCH ${CMAKE_MATCH_3})
      set(VERSION_INFO ${CMAKE_MATCH_4})
    else()
      # default values
      set(VERSION_MAJOR "0")
      set(VERSION_MINOR "0")
      set(VERSION_PATCH "0")
    endif()
  endif()

  set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endmacro(get_version_from_git)


macro(print_debug_messages)
  message(DEBUG "CMAKE_C_FLAGS is: ${CMAKE_C_FLAGS}")
  message(DEBUG "CMAKE_C_FLAGS_DEBUG is: ${CMAKE_C_FLAGS_DEBUG}")
  message(DEBUG "CMAKE_C_FLAGS_RELEASE is: ${CMAKE_C_FLAGS_RELEASE}")
  message(DEBUG "CMAKE_C_FLAGS_RELWITHDEBINFO is: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  message(DEBUG "CMAKE_C_FLAGS_MINSIZEREL is: ${CMAKE_C_FLAGS_MINSIZEREL}")
  message(DEBUG "CMAKE_CXX_FLAGS is: ${CMAKE_CXX_FLAGS}")
  message(DEBUG "CMAKE_CXX_FLAGS_DEBUG is: ${CMAKE_CXX_FLAGS_DEBUG}")
  message(DEBUG "CMAKE_CXX_FLAGS_RELEASE is: ${CMAKE_CXX_FLAGS_RELEASE}")
  message(DEBUG "CMAKE_CXX_FLAGS_RELWITHDEBINFO is: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  message(DEBUG "CMAKE_CXX_FLAGS_MINSIZEREL is: ${CMAKE_CXX_FLAGS_MINSIZEREL}")
  message(DEBUG "Build type: ${CMAKE_BUILD_TYPE}")
  message(DEBUG "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE})
  message(DEBUG "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS})
endmacro(print_debug_messages)
