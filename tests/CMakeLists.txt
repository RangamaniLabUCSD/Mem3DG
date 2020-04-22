cmake_minimum_required(VERSION 3.14)

include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY    https://github.com/google/googletest.git
    GIT_TAG           release-1.10.0
    SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
)
# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${googletest_SOURCE_DIR}
                     ${googletest_BINARY_DIR}
                     EXCLUDE_FROM_ALL)
endif()


# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

set_property(TARGET gtest PROPERTY CXX_STANDARD 11)

# Build the tests
set(TEST_SRCS
  src/main_test.cpp
  src/mapping_test.cpp
)

add_executable(ddgsolver-test "${TEST_SRCS}")
target_link_libraries(ddgsolver-test ddg_solver gtest_main)

add_test(NAME ddgtests COMMAND ddgsolver-test)