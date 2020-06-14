include(FetchContent)

add_subdirectory(pcg)

# ##############################################################################
# GET GEOMETRY CENTRAL
# ##############################################################################
FetchContent_Declare(
  geometrycentrallib
  GIT_REPOSITORY https://github.com/ctlee/geometry-central.git
  GIT_TAG origin/data_access
  GIT_SHALLOW TRUE
  SOURCE_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-src"
  BINARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-build"
)
FetchContent_GetProperties(geometrycentrallib)
if(NOT geometrycentrallib_POPULATED)
  FetchContent_Populate(geometrycentrallib)
  add_subdirectory(
    ${geometrycentrallib_SOURCE_DIR} ${geometrycentrallib_BINARY_DIR}
  )
  set_target_properties(geometry-central PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
endif()

# ##############################################################################
# GET polyscope
# ##############################################################################
FetchContent_Declare(
  polyscopelib
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_SHALLOW TRUE
  SOURCE_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/polyscope-src"
  BINARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/polyscope-build"
)
FetchContent_GetProperties(polyscopelib)
if(NOT polyscopelib_POPULATED)
  FetchContent_Populate(polyscopelib)
  add_subdirectory(${polyscopelib_SOURCE_DIR} ${polyscopelib_BINARY_DIR})
endif()
