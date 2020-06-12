cmake_minimum_required(VERSION 3.11)

if(MSVC)
    set(CMAKE_C_FLAGS_DEBUG_INIT "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_C_FLAGS_MINSIZEREL_INIT     "/MT /O1 /Ob1 /D NDEBUG")
    set(CMAKE_C_FLAGS_RELEASE_INIT        "/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "/MT /Zi /O2 /Ob1 /D NDEBUG")
else()
    set(CMAKE_C_FLAGS_RELEASE_INIT        "-O2 -DNDEBUG")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/MT /O1 /Ob1 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "/MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/MT /Zi /O2 /Ob1 /D NDEBUG")
else()
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "-O2 -DNDEBUG")
endif()

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

option(BUILD_PYDDG "Build the python extensions?" ON)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/velocity_verlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bending_energy.cpp
)

add_library(ddg_solver SHARED ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(ddg_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ddg_solver PUBLIC geometry-central pcg::pcg polyscope)


##########################
# Main DDG target
##########################
add_executable(ddg tests/main.cpp)
target_link_libraries(ddg PRIVATE ddg_solver)

if(BUILD_PYDDG)
    add_subdirectory(pyddglib)
endif()
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
