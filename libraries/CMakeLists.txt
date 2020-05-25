
include(FetchContent)

add_subdirectory(pcg)

##########################
# GET PYBIND11
##########################
if(BUILD_PYDDG)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY  https://github.com/pybind/pybind11.git
        GIT_TAG         v2.5.0
        GIT_SHALLOW     TRUE
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
    )
    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
      FetchContent_Populate(pybind11)
      add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
endif()

##########################
# GET GEOMETRY CENTRAL
##########################
FetchContent_Declare(
    geometrycentrallib
    GIT_REPOSITORY  https://github.com/ctlee/geometry-central.git
    GIT_TAG         origin/data_access
    GIT_SHALLOW     TRUE
    SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-src"
    BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-build"
)
FetchContent_GetProperties(geometrycentrallib)
if(NOT geometrycentrallib_POPULATED)
  FetchContent_Populate(geometrycentrallib)
  add_subdirectory(${geometrycentrallib_SOURCE_DIR} ${geometrycentrallib_BINARY_DIR})
endif()
