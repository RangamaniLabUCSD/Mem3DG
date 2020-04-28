cmake_minimum_required(VERSION 3.8.0)

project(pcg VERSION 0.98.1 LANGUAGES CXX)

include(GNUInstallDirs)

add_library(pcg INTERFACE)
add_library(pcg::pcg ALIAS pcg)

target_compile_features(pcg INTERFACE cxx_std_11)

target_include_directories(pcg
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)