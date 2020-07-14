# ##############################################################################
# GET PYBIND11
# ##############################################################################
if(BUILD_PYDDG)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.5.0
    GIT_SHALLOW TRUE
    SOURCE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/pybind11-src"
    BINARY_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/pybind11-build"
  )
  FetchContent_GetProperties(pybind11)
  if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
  endif()
endif()

# ddgsolver static library
add_library(
  ddgsolver_static STATIC EXCLUDE_FROM_ALL $<TARGET_OBJECTS:ddgsolver_objlib>
)
target_link_libraries(ddgsolver_static PUBLIC ddgsolver_objlib)
set_target_properties(ddgsolver_static PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden)

#set(PYDDG_SOURCES "src/pyddg.cpp")
set(PYDDG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pyddg.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ddgsolver.cpp
)

pybind11_add_module(pyddg MODULE ${PYDDG_SOURCES})
target_link_libraries(pyddg PRIVATE ddgsolver_static polyscope)
