

# WORKING PKG CONFIG!
# find_package(PkgConfig QUIET)
# if(PkgConfig_FOUND)
#   # message(DEBUG "PKGCONFIG FOUND")
#   pkg_check_modules(_netcdf-cxx4 QUIET netcdf-cxx4 IMPORTED_TARGET)
#   if(_netcdf-cxx4_FOUND)
#     # message(DEBUG "netcdf-cxx4 found by pkgconfig")
#     # Forward the variables in a consistent way.
#     set(netcdf-cxx4_FOUND "${_netcdf-cxx4_FOUND}")
#     set(netcdf-cxx4_INCLUDE_DIRS "${_netcdf-cxx4_INCLUDE_DIRS}")
#     set(netcdf-cxx4_LIBRARIES "${_netcdf-cxx4_LIBRARIES}")
#     set(netcdf-cxx4_VERSION "${_netcdf-cxx4_VERSION}")
#     if(NOT TARGET NetCDF::NetCDF-cxx4)
#       add_library(NetCDF::NetCDF-cxx4 INTERFACE IMPORTED)
#       set_target_properties(
#         NetCDF::NetCDF-cxx4 PROPERTIES INTERFACE_LINK_LIBRARIES
#                                        "PkgConfig::_netcdf-cxx4"
#       )
#     endif()
#     # Skip the rest of the logic in this file.
#     return()
#   endif()
# endif()

find_path(
  netcdf-cxx4_INCLUDE_DIR
  NAMES netcdf
  DOC "netcdf-cxx4 include directories"
)
mark_as_advanced(netcdf-cxx4_INCLUDE_DIR)

find_library(
  netcdf-cxx4_LIBRARY_RELEASE
  NAMES netcdf-cxx4 netcdf_c++4
  DOC "netcdf-cxx4 release library"
)

find_library(
  netcdf-cxx4_LIBRARY_DEBUG
  NAMES netcdf-cxx4 netcdf_c++4
  PATH_SUFFIXES debug/lib
  DOC "netcdf-cxx4 debug library"
)
include(SelectLibraryConfigurations)
select_library_configurations(netcdf-cxx4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  netcdf-cxx4 REQUIRED_VARS netcdf-cxx4_LIBRARY netcdf-cxx4_INCLUDE_DIR
  # VERSION_VAR NetCDF_VERSION
)

if(netcdf-cxx4_FOUND)
  set(netcdf-cxx4_INCLUDE_DIRS "${netcdf-cxx4_INCLUDE_DIR}")
  set(netcdf-cxx4_LIBRARIES "${netcdf-cxx4_LIBRARY}")

  if(NOT TARGET NetCDF::NetCDF-cxx4)
    add_library(NetCDF::NetCDF-cxx4 UNKNOWN IMPORTED)
    set_target_properties(
      NetCDF::NetCDF-cxx4
      PROPERTIES IMPORTED_LOCATION "${netcdf-cxx4_LIBRARY}"
                 INTERFACE_INCLUDE_DIRECTORIES "${netcdf-cxx4_INCLUDE_DIR}"
    )
  endif()
endif()
