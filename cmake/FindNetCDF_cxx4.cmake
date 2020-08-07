
# From VTK
# https://gitlab.kitware.com/vtk/vtk/blob/master/CMake/FindNetCDF.cmake

#[==[
Provides the following variables:

  * `NetCDF_cxx4_FOUND`: Whether NetCDF was found or not.
  * `NetCDF_cxx4_INCLUDE_DIRS`: Include directories necessary to use NetCDF.
  * `NetCDF_cxx4_LIBRARIES`: Libraries necessary to use NetCDF.
  * `NetCDF_VERSION`: The version of NetCDF found.
  * `NetCDF::NetCDF_cxx4`: A target to use with `target_link_libraries`.
#]==]

# Try to find a CMake-built NetCDF.
find_package(netCDF_cxx4 CONFIG QUIET)
if (netCDF_cxx4_FOUND)
  # Forward the variables in a consistent way.
  set(NetCDF_cxx4_FOUND "${netCDF_cxx4_FOUND}")
  set(NetCDF_cxx4_INCLUDE_DIRS "${netCDF_cxx4_INCLUDE_DIR}")
  set(NetCDF_cxx4_LIBRARIES "${netCDF_cxx4_LIBRARIES}")
  set(NetCDF_cxx4_VERSION "${NetCDFVersion}")
  if (NOT TARGET NetCDF::NetCDF_cxx4)
    add_library(NetCDF::NetCDF_cxx4 INTERFACE IMPORTED)
    if (TARGET "netCDF::netcdf_cxx4")
      # 4.7.3
      set_target_properties(NetCDF::NetCDF_cxx4 PROPERTIES
        INTERFACE_LINK_LIBRARIES "netCDF::netcdf_cxx4")
    elseif (TARGET "netcdf_cxx4")
      set_target_properties(NetCDF::NetCDF_cxx4 PROPERTIES
        INTERFACE_LINK_LIBRARIES "netcdf_cxx4")
    else ()
      set_target_properties(NetCDF::NetCDF_cxx4 PROPERTIES
        INTERFACE_LINK_LIBRARIES "${netCDF_cxx4_LIBRARIES}")
    endif ()
  endif ()
  # Skip the rest of the logic in this file.
  return ()
endif ()

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(_NetCDF QUIET netcdf_cxx4 IMPORTED_TARGET)
  if (_NetCDF_FOUND)
    # Forward the variables in a consistent way.
    set(NetCDF_cxx4_FOUND "${_NetCDF_FOUND}")
    set(NetCDF_cxx4_INCLUDE_DIRS "${_NetCDF_INCLUDE_DIRS}")
    set(NetCDF_cxx4_LIBRARIES "${_NetCDF_LIBRARIES}")
    set(NetCDF_cxx4_VERSION "${_NetCDF_VERSION}")
    if (NOT TARGET NetCDF::NetCDF_cxx4)
      add_library(NetCDF::NetCDF_cxx4 INTERFACE IMPORTED)
      set_target_properties(NetCDF::NetCDF_cxx4 PROPERTIES
        INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDF")
    endif ()
    # Skip the rest of the logic in this file.
    return ()
  endif ()
endif ()

find_path(NetCDF_cxx4_INCLUDE_DIR
  NAMES netcdf
  DOC "netcdf-cxx4 include directories")
mark_as_advanced(NetCDF_cxx4_INCLUDE_DIR)

find_library(NetCDF_cxx4_LIBRARY
  NAMES netcdf_c++4
  DOC "netcdf-cxx4 library")
mark_as_advanced(NetCDF_cxx4_LIBRARY)

if (NetCDF_cxx4_INCLUDE_DIR)
  file(STRINGS "${NetCDF_cxx4_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
    REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
  string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
  set(NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
  unset(_netcdf_version_major)
  unset(_netcdf_version_minor)
  unset(_netcdf_version_patch)
  unset(_netcdf_version_note)
  unset(_netcdf_version_lines)
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF_cxx4
  REQUIRED_VARS NetCDF_cxx4_LIBRARY NetCDF_cxx4_INCLUDE_DIR
  VERSION_VAR NetCDF_VERSION)

if (NetCDF_cxx4_FOUND)
  set(NetCDF_cxx4_INCLUDE_DIRS "${NetCDF_cxx4_INCLUDE_DIR}")
  set(NetCDF_cxx4_LIBRARIES "${NetCDF_cxx4_LIBRARY}")

  if (NOT TARGET NetCDF::NetCDF_cxx4)
    add_library(NetCDF::NetCDF_cxx4 UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_cxx4 PROPERTIES
      IMPORTED_LOCATION "${NetCDF_cxx4_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_cxx4_INCLUDE_DIR}")
  endif ()
endif ()
