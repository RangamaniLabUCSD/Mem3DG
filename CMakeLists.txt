cmake_minimum_required(VERSION 3.11)

# Disable in source builds
# Note that this is an undocumented option which may not do anything
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(DDG VERSION 0.0.1 LANGUAGES CXX)

# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize output directory locations
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Required to tell MSVC to export symbols
if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Check if DDG is being used directly or via add_subdirectory
set(DDG_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(DDG_MASTER_PROJECT ON)
endif()

# find_package(Git REQUIRED)
include(FetchContent)

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

##########################
# BUNDLED LIBRARIES
##########################
add_subdirectory(libraries)

##########################
# DDG SOLVER LIBRARY
##########################
add_subdirectory(include) # Get file sources

set(DDG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/icosphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bending_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stretching_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pressure_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/damping_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stochastic_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stormer_verlet.cpp
)

add_library(ddg_solver SHARED ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(ddg_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ddg_solver PUBLIC geometry-central pcg::pcg)


##########################
# Main DDG target
##########################
add_executable(ddg tests/main.cpp)
target_link_libraries(ddg PRIVATE ddg_solver)

##########################
# GET polyscope
##########################
FetchContent_Declare(
    polyscopelib
    GIT_REPOSITORY  https://github.com/nmwsharp/polyscope.git
    GIT_SHALLOW     TRUE
    SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/polyscope-src"
    BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/polyscope-build"
)
FetchContent_GetProperties(polyscopelib)
if(NOT polyscopelib_POPULATED)
  FetchContent_Populate(polyscopelib)
  add_subdirectory(${polyscopelib_SOURCE_DIR} ${polyscopelib_BINARY_DIR})
endif()

target_link_libraries(ddg PRIVATE polyscope)

enable_testing()
add_subdirectory(tests)
