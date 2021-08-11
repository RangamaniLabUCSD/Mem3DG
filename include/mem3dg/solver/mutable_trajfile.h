// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

/**
 * @file  mutable_MutableTrajFile.h
 * @brief Netcdf trajectory output support
 *
 */

#pragma once

#ifdef MEM3DG_WITH_NETCDF

#include <exception>
#include <iostream>
#include <netcdf>
#include <vector>

#include <Eigen/Core>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "mem3dg/macros.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/trajfile_constants.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;

/**
 * @class MutableTrajFile
 * @brief Trajectory interface to help with manipulating trajectories
 *
 */
class DLL_PUBLIC MutableTrajFile {
private:
  /* Do not allow definition of NcFile involving copying any NcFile or
    NcGroup. Because the destructor closes the file and releases al resources
    such an action could leave NcFile objects in an invalid state from NcFile.h
  */
  // MutableTrajFile *operator=(const MutableTrajFile &rhs) = delete;
  MutableTrajFile(const MutableTrajFile &rhs) = delete;

public:
  using NcFile = nc::NcFile;
  using NcException = nc::exceptions::NcException;

#pragma region named_constructors
  /**
   * @brief Open a new file and populate it with the convention
   *
   * @param filename  Filename to save to
   * @param replace   Whether to replace an existing file or exit
   *
   * @exception netCDF::exceptions::NcExist File already exists and
   * replace/overwrite flag is not specified.
   *
   * @return MutableTrajFile helper object to manipulate the bound NetCDF file.
   */
  static MutableTrajFile newFile(const std::string &filename,
                                 bool replace = false) {
    if (replace)
      return MutableTrajFile(filename, NcFile::replace);
    else
      return MutableTrajFile(filename, NcFile::newFile);
  };

  /**
   * @brief Open a new file and populate it with the convention
   *
   * @param filename  Filename to save to
   * @param replace   Whether to replace an existing file or exit
   *
   * @exception netCDF::exceptions::NcExist File already exists and
   * replace/overwrite flag is not specified.
   *
   * @return MutableTrajFile helper object to manipulate the bound NetCDF file.
   */
  static MutableTrajFile newDisklessFile(const std::string &filename) {
    return MutableTrajFile(filename, NC_NETCDF4 | NC_CLOBBER | NC_DISKLESS);
  };

  /**
   * @brief Open an existing NetCDF file in read/write mode
   *
   * @param filename  Filename of interest
   *
   * @exception std::runtime_error if file does not conform to the convention
   * @exception netCDF::exceptions::* If file does not exist
   *
   * @return MutableTrajFile helper object to manipulate the bound NetCDF file.
   */
  static MutableTrajFile openRW(const std::string &filename) {
    return MutableTrajFile(filename, NcFile::write);
  }

  /**
   * @brief Open an existing NetCDF file in read only mode
   *
   * @param filename  Filename of interest
   *
   * @exception std::runtime_error if file does not conform to the convention
   * @exception netCDF::exceptions::* If file does not exist
   *
   * @return MutableTrajFile helper object to manipulate the bound NetCDF file.
   */
  static MutableTrajFile openReadOnly(const std::string &filename) {
    return MutableTrajFile(filename, NcFile::read);
  };
#pragma endregion named_constructors

  /// Default constructor
  MutableTrajFile() : writeable(false), fd(nullptr){};

  /// Default copy constructor
  MutableTrajFile(MutableTrajFile &&rhs) = default;

  /**
   * @brief Destructor frees the bound NcFile
   */
  ~MutableTrajFile() { delete fd; };

#pragma region initialization_helpers
  /**
   * @brief Open an existing trajectory file for reading/writing
   *
   * @param filename  Path to file to open
   * @param fMode     Mode to open file with
   */
  void open(const std::string &filename, const NcFile::FileMode fMode) {
    if ((fd != nullptr) && (fMode != NcFile::read)) {
      mem3dg_runtime_error("Cannot open an already opened ...");
    }

    fd = new NcFile(filename, fMode);
    writeable = fMode != NcFile::read;
    check_metadata();

    frame_dim = fd->getDim(FRAME_NAME);
  }

  /**
   * @brief Create a new file with NcFile file modes (from NetCDF-C++4)
   *
   * @param filename  Path to file to create
   * @param fMode     Mode to create the file
   */
  void createNewFile(const std::string &filename,
                     const NcFile::FileMode fMode) {
    if (fd != nullptr) {
      mem3dg_runtime_error("Cannot open an already open ...");
    }

    writeable = true;

    fd = new NcFile(filename, fMode);
    initializeConventions();
  }

  /**
   * @brief Create a New File object with NetCDF-C modes
   *
   * @param filename    Path to file to create
   * @param ncFileMode  Mode to create the file
   */
  void createNewFile(const std::string &filename, const int ncFileMode) {
    if (fd != nullptr) {
      mem3dg_runtime_error("Cannot open an already open ...");
    }

    writeable = true;

    fd = new NcFile();
    fd->create(filename, ncFileMode);
    initializeConventions();
  }
#pragma endregion initialization_helpers

  /**
   * @brief Netcdf file frame reader
   *
   * @param frame reference to the frame index
   *
   */
  void getNcFrame(int &frame) const {
    int maxFrame = getNextFrameIndex() - 1;
    if (frame > maxFrame || frame < -(maxFrame + 1)) {
      mem3dg_runtime_error("Snapshot frame exceed limiting frame index!");
    } else if (frame < 0) {
      frame = frame + maxFrame + 1;
    }
  }

  /**
   * @brief Check if the bound NcFile was opened in write mode
   *
   * @return True if writable
   */
  inline bool isWriteable() { return writeable; };

  inline std::size_t getNextFrameIndex() const { return frame_dim.getSize(); };

  void writeTopoFrame(const std::size_t idx, const EigenVectorX3ur &data);

  EigenVectorX3ur readTopoFrame(const std::size_t idx);

  void writeTime(const std::size_t idx, const double time);

  /**
   * @brief Validate whether or not the metadata follows convention
   *
   * @return True if okay. False otherwise.
   */
  bool check_metadata();

private:
  /**
   * @brief Private constructor for opening or creating a new NetCDF file.
   *
   * The metadata is checked for consistency with the convention.
   *
   * @param filename Path to file of interest
   * @param fMode    Mode to open/create file with
   */
  MutableTrajFile(const std::string &filename, const NcFile::FileMode fMode)
      : filename(filename), fd(nullptr), writeable(fMode != NcFile::read) {
    if (fMode == NcFile::read || fMode == NcFile::write)
      open(filename, fMode);
    else
      createNewFile(filename, fMode);
  }

  /**
   * @brief Private constructor for new file using ncflags directly
   *
   * Note that this only creates new files!
   *
   * @param filename      Path to file of interest
   * @param ncFileFlags   Mode to create file with
   */
  MutableTrajFile(const std::string &filename, const int ncFileFlags)
      : filename(filename), fd(nullptr), writeable(true) {
    createNewFile(filename, ncFileFlags);
  }

  /**
   * @brief Initialize a new file with the given conventions
   */
  void initializeConventions() {
    // initialize data
    fd->putAtt(CONVENTIONS_NAME, CONVENTIONS_VALUE);
    fd->putAtt(CONVENTIONS_VERSION_NAME, CONVENTIONS_VERSION_VALUE);

    frame_dim = fd->addDim(FRAME_NAME);

    time_var = fd->addVar(TIME_VAR, netCDF::ncDouble, {frame_dim});
    time_var.putAtt(UNITS, TIME_UNITS);

    uint_array_t = fd->addVlenType(UINT_ARR, nc::ncUint);
    double_array_t = fd->addVlenType(DOUBLE_ARR, nc::ncDouble);

    topo_var = fd->addVar(TOPO_VAR, uint_array_t, {frame_dim});
  }

  /// Bound NcFile
  NcFile *fd;

  // Save dimensions
  nc::NcDim frame_dim;

  /// Variable length type for topology
  nc::NcVlenType uint_array_t;
  nc::NcVlenType double_array_t;

  /// Variable for storing time
  nc::NcVar time_var;
  /// Vlen variable for topology
  nc::NcVar topo_var;

  /// Filepath to file
  std::string filename;
  /// Writeable status
  bool writeable;
};
} // namespace solver
} // namespace mem3dg
#endif
