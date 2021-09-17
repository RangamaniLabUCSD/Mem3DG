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

  ///////////////////////////////////////////////
  // Reintroduce when Netcdf 4.3.1 is standard
  // Ubuntu 18.04 uses 4.3.0
  ///////////////////////////////////////////////
  //   /**
  //    * @brief Open a new file and populate it with the convention
  //    *
  //    * @param filename  Filename to save to
  //    * @param replace   Whether to replace an existing file or exit
  //    *
  //    * @exception netCDF::exceptions::NcExist File already exists and
  //    * replace/overwrite flag is not specified.
  //    *
  //    * @return MutableTrajFile helper object to manipulate the bound NetCDF
  //    file.
  // ;   */
  //   static MutableTrajFile newDisklessFile(const std::string &filename) {
  //     return MutableTrajFile(filename, NC_NETCDF4 | NC_CLOBBER |
  //     NC_DISKLESS);
  //   }

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
  ~MutableTrajFile() {
    if (fd != nullptr) {
      fd->close();
    }
    delete fd;
  };

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

    parameter_group = fd->getGroup(PARAM_GROUP_NAME);
    traj_group = fd->getGroup(TRAJ_GROUP_NAME);

    frame_dim = traj_group.getDim(FRAME_NAME);

    // uint_array_t = traj_group.getType(UINT_ARR);
    // double_array_t = traj_group.getType(DOUBLE_ARR);

    time_var = traj_group.getVar(TIME_VAR);
    topo_var = traj_group.getVar(TOPO_VAR);
    coord_var = traj_group.getVar(COORD_VAR);
    phi_var = traj_group.getVar(PHI_VAR);
    vel_var = traj_group.getVar(VEL_VAR);
    extF_var = traj_group.getVar(EXTF_VAR);
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

  ///////////////////////////////////////////////
  // Reintroduce when Netcdf 4.3.1 is standard
  // Ubuntu 18.04 uses 4.3.0
  ///////////////////////////////////////////////
  // /**
  //  * @brief Create a New File object with NetCDF-C modes
  //  *
  //  * @param filename    Path to file to create
  //  * @param ncFileMode  Mode to create the file
  //  */
  // void createNewFile(const std::string &filename, const int ncFileMode) {
  //   if (fd != nullptr) {
  //     mem3dg_runtime_error("Cannot open an already opened file.");
  //   }

  //   writeable = true;

  //   fd = new NcFile();
  //   fd->create(filename, ncFileMode);
  //   initializeConventions();
  // }

#pragma endregion initialization_helpers

  void sync() { fd->sync(); }

  /**
   * @brief Close
   *
   */
  void close() {
    if (fd == nullptr) {
      mem3dg_runtime_error("Cannot close an unopened trajectory file.");
    }
    fd->sync();
    fd->close();

    // Reset object state
    fd = nullptr;
    writeable = false;

    traj_group = nc::NcGroup{};
    parameter_group = nc::NcGroup{};

    frame_dim = nc::NcDim{};
    uint_array_t = nc::NcVlenType{};
    double_array_t = nc::NcVlenType{};
    time_var = nc::NcVar{};
    topo_var = nc::NcVar{};
    coord_var = nc::NcVar{};
    phi_var = nc::NcVar{};
    vel_var = nc::NcVar{};
    extF_var = nc::NcVar{};
    filename = "";
  }

  /**
   * @brief Check if the bound NcFile was opened in write mode
   *
   * @return True if writable
   */
  bool isWriteable() { return writeable; };

  /**
   * @brief Get the max number of frames in the trajectory
   *
   * @return std::size_t Total number of frames
   */
  std::size_t nFrames() const { return frame_dim.getSize(); };

  /**
   * @brief Netcdf file frame reader
   *
   * @param frame reference to the frame index
   *
   */
  void getNcFrame(int &frame) const {
    int maxFrame = nFrames() - 1;
    if (frame > maxFrame || frame < -(maxFrame + 1)) {
      mem3dg_runtime_error("Snapshot frame exceed limiting frame index!");
    } else if (frame < 0) {
      frame = frame + maxFrame + 1;
    }
  }

#pragma region read_write
  /**
   * @brief Write the topology for a frame
   *
   * @param idx   Index of the frame
   * @param data  Topology matrix
   */
  void writeTopology(const std::size_t idx, const EigenVectorX3ur &data) {
    writeVar<std::uint32_t, 3>(topo_var, idx, data);
  }

  /**
   * @brief Write the topology for a frame
   *
   * @param idx   Index of the frame
   * @param data  Surface mesh
   */
  void writeTopology(const std::size_t idx, gc::SurfaceMesh &mesh) {
    writeVar<std::uint32_t, 3>(
        topo_var, idx,
        EigenVectorX3ur{mesh.getFaceVertexMatrix<std::uint32_t>()});
  }

  /**
   * @brief Get the Topology object
   *
   * @param idx
   * @return EigenVectorX3ur
   */
  EigenVectorX3ur getTopology(const std::size_t idx) const {
    return getVar<std::uint32_t, POLYGON_ORDER>(topo_var, idx);
  }

  /**
   * @brief Write the coordinates for a frame
   *
   * @param idx   Index of the frame
   * @param data  Coordinate matrix
   */
  void writeCoords(const std::size_t idx, const EigenVectorX3dr &data) {
    writeVar<double, 3>(coord_var, idx, data);
  }

  /**
   * @brief Write the coordinates for a frame
   *
   * @param idx   Index of the frame
   * @param data  Vertex position geometry
   */
  void writeCoords(const std::size_t idx,
                   const gc::VertexPositionGeometry &data) {
    writeVar<gc::Vertex>(coord_var, idx, data.inputVertexPositions);
  }

  /**
   * @brief Get the coordinates of a given frame
   *
   * @param idx               Index of the frame
   * @return EigenVectorX3dr  Coordinates data
   */
  EigenVectorX3dr getCoords(const std::size_t idx) {
    return getVar<double, SPATIAL_DIMS>(coord_var, idx);
  }

  /**
   * @brief Write the protein density for a frame
   *
   * @param idx   Index of the frame
   * @param data  Vertex position geometry
   */
  void writeProteinDensity(const std::size_t idx,
                           const gc::MeshData<gc::Vertex, double> &data) {
    writeVar(phi_var, idx, data);
  }

  /**
   * @brief Get the protein density of a given frame
   *
   * @param idx               Index of the frame
   * @return EigenVectorX3dr  Coordinates data
   */
  EigenVectorX1d getProteinDensity(const std::size_t idx) {
    return getVar1d<double>(phi_var, idx);
  }

  /**
   * @brief Write the velocities for a frame
   *
   * @param idx   Index of the frame
   * @param data  Velocity matrix
   */
  void writeVelocity(const std::size_t idx, const EigenVectorX3dr &data) {
    writeVar<double, 3>(vel_var, idx, data);
  }

  /**
   * @brief Write the velocities for a frame
   *
   * @param idx   Index of the frame
   * @param data  Vertex velocities
   */
  void writeVelocity(const std::size_t idx,
                     const gcs::VertexData<gc::Vector3> &data) {
    writeVar<gc::Vertex>(vel_var, idx, data);
  }

  /**
   * @brief Get the velocities of a given frame
   *
   * @param idx               Index of the frame
   * @return EigenVectorX3dr  Velocity data
   */
  EigenVectorX3dr getVelocity(const std::size_t idx) const {
    return getVar<double, SPATIAL_DIMS>(vel_var, idx);
  }

  /**
   * @brief Write the external force field for a frame
   *
   * @param idx   Index of the frame
   * @param data  Velocity matrix
   */
  void writeExternalForce(const std::size_t idx, const EigenVectorX3dr &data) {
    writeVar<double, 3>(extF_var, idx, data);
  }

  /**
   * @brief Write the external force field for a frame
   *
   * @param idx   Index of the frame
   * @param data  Vertex velocities
   */
  void writeExternalForce(const std::size_t idx,
                          const gcs::VertexData<gc::Vector3> &data) {
    writeVar<gc::Vertex>(extF_var, idx, data);
  }

  /**
   * @brief Get the external force field of a given frame
   *
   * @param idx               Index of the frame
   * @return EigenVectorX3dr  Velocity data
   */
  EigenVectorX3dr getExternalForce(const std::size_t idx) const {
    return getVar<double, SPATIAL_DIMS>(extF_var, idx);
  }

  /**
   * @brief Write the time of the trajectory
   *
   * @param idx     Index of the frame
   * @param time    Time
   */
  void writeTime(const std::size_t idx, const double time) {
    writeVar(time_var, idx, time);
  }

  /**
   * @brief Get the time
   *
   * @param idx      Index
   * @return double  Time
   */
  double getTime(const std::size_t idx) const {
    return getVar<double>(time_var, idx);
  }

#pragma endregion read_write
  /**
   * @brief Validate whether or not the metadata follows convention
   *
   * @return True if okay. False otherwise.
   */
  bool check_metadata();

private:
  // template <typename SCALAR, std::size_t k, typename T>
  // void writeVar(nc::NcVar &var, const std::size_t idx, const T &data,
  //               std::function<EigenVectorXkr_T<SCALAR, k>(T)>
  //               &&toEigenVector) {
  //   writeVar(var, idx, toEigenVector(data));
  // }

  // /**
  //  * @brief Write double MeshData to variable
  //  *
  //  * @tparam E      Typename of the Mesh Element
  //  * @param var     Variable to write to
  //  * @param idx     Index
  //  * @param data    Data
  //  */
  // template <typename E>
  // void writeVar(nc::NcVar &var, const std::size_t idx,
  //               const gc::MeshData<E, double> &data) {
  //   writeVar<E, double>(var, idx, data);
  // }

  /**
   * @brief Write gc::Vector3 MeshData to variable
   *
   * @tparam E      Typename of the Mesh Element
   * @param var     Variable to write to
   * @param idx     Index
   * @param data    Data
   */
  template <typename E>
  void writeVar(nc::NcVar &var, const std::size_t idx,
                const gc::MeshData<E, gc::Vector3> &data) {
    writeVar<double, 3>(var, idx, EigenMap<double, 3>(data));
  }

  /**
   * @brief Write MeshData storing a primitive type to a variable
   *
   * @tparam E      Typename of the Mesh Element
   * @tparam T      Typename of the data type
   * @param var     Variable to write to
   * @param idx     Index
   * @param data    Data
   */
  template <typename E, typename T,
            typename = std::enable_if_t<std::is_fundamental<T>::value>>
  void writeVar(nc::NcVar &var, const std::size_t idx,
                const gc::MeshData<E, T> &data) {
    writeVar(var, idx, data.raw());
  }

  /**
   * @brief
   *
   * @tparam T
   * @tparam k
   * @param var
   * @param idx
   * @param data
   */
  template <typename T, int k>
  void writeVar(nc::NcVar &var, const std::size_t idx,
                const EigenVectorXkr_T<T, k> &data) {
    if (!writeable)
      mem3dg_runtime_error("Cannot write to read only file.");

    nc_vlen_t vlenData;
    vlenData.len = data.size();
    vlenData.p = const_cast<T *>(data.data());

    var.putVar({idx}, &vlenData);
  }

  /**
   * @brief
   *
   * @tparam T
   * @param var
   * @param idx
   * @param data
   */
  template <typename T>
  void writeVar(nc::NcVar &var, const std::size_t idx,
                const EigenVectorX1_T<T> &data) {
    if (!writeable)
      mem3dg_runtime_error("Cannot write to read only file.");

    nc_vlen_t vlenData;
    vlenData.len = data.size();
    vlenData.p = const_cast<T *>(data.data());

    var.putVar({idx}, &vlenData);
  }

  template <typename T,
            typename = std::enable_if_t<std::is_fundamental<T>::value>>
  void writeVar(nc::NcVar &var, const std::size_t idx, const T data) {
    if (!writeable)
      mem3dg_runtime_error("Cannot write to read only file.");
    var.putVar({idx}, &data);
  }

  template <typename T, std::size_t k>
  EigenVectorXkr_T<T, k> getVar(const nc::NcVar &var,
                                const std::size_t idx) const {
    assert(idx < nFrames());

    nc_vlen_t vlenData;
    var.getVar({idx}, &vlenData);

    // Initialize an Eigen object and copy the data over
    EigenVectorXkr_T<T, k> vec(vlenData.len / k, k);
    // Bind to nc_vlen_t memory and copy data over
    vec = AlignedEigenMap_T<T, k, Eigen::RowMajor>(static_cast<T *>(vlenData.p),
                                                   vlenData.len / k, k);
    return vec;
  }

  template <typename T>
  EigenVectorX1_T<T> getVar1d(const nc::NcVar &var,
                                const std::size_t idx) const {
    assert(idx < nFrames());

    nc_vlen_t vlenData;
    var.getVar({idx}, &vlenData);

    // Initialize an Eigen object and copy the data over
    EigenVectorX1_T<T> vec(vlenData.len);
    // Bind to nc_vlen_t memory and copy data over
    vec = AlignedEigenMap_T<T, 1, Eigen::ColMajor>(static_cast<T *>(vlenData.p), vlenData.len);
    return vec;
  }


  template <typename T,
            typename = std::enable_if_t<std::is_fundamental<T>::value>>
  T getVar(const nc::NcVar &var, const std::size_t idx) const {
    assert(idx < nFrames());

    T data;
    var.getVar({idx}, &data);

    return data;
  }

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

  ///////////////////////////////////////////////
  // Reintroduce when Netcdf 4.3.1 is standard
  // Ubuntu 18.04 uses 4.3.0
  ///////////////////////////////////////////////
  // /**
  //  * @brief Private constructor for new file using ncflags directly
  //  *
  //  * Note that this only creates new files!
  //  *
  //  * @param filename      Path to file of interest
  //  * @param ncFileFlags   Mode to create file with
  //  */
  // MutableTrajFile(const std::string &filename, const int ncFileFlags)
  //     : filename(filename), fd(nullptr), writeable(true) {
  //   createNewFile(filename, ncFileFlags);
  // }

  /**
   * @brief Initialize a new file with the given conventions
   */
  void initializeConventions() {

    const int compression_level = 5;

    // initialize data
    fd->putAtt(CONVENTIONS_NAME, CONVENTIONS_VALUE);
    fd->putAtt(CONVENTIONS_VERSION_NAME, CONVENTIONS_VERSION_VALUE);

    parameter_group = fd->addGroup(PARAM_GROUP_NAME);
    traj_group = fd->addGroup(TRAJ_GROUP_NAME);

    frame_dim = traj_group.addDim(FRAME_NAME);

    time_var = traj_group.addVar(TIME_VAR, netCDF::ncDouble, {frame_dim});
    time_var.putAtt(UNITS, TIME_UNITS);
    time_var.setCompression(true, true, compression_level);

    uint_array_t = traj_group.addVlenType(UINT_ARR, nc::ncUint);
    double_array_t = traj_group.addVlenType(DOUBLE_ARR, nc::ncDouble);

    topo_var = traj_group.addVar(TOPO_VAR, uint_array_t, {frame_dim});
    topo_var.setCompression(true, true, compression_level);
    coord_var = traj_group.addVar(COORD_VAR, double_array_t, {frame_dim});
    coord_var.setCompression(true, true, compression_level);
    phi_var = traj_group.addVar(PHI_VAR, double_array_t, {frame_dim});
    phi_var.setCompression(true, true, compression_level);
    vel_var = traj_group.addVar(VEL_VAR, double_array_t, {frame_dim});
    vel_var.setCompression(true, true, compression_level);
    extF_var = traj_group.addVar(EXTF_VAR, double_array_t, {frame_dim});
    extF_var.setCompression(true, true, compression_level);
  }

  /// Bound NcFile
  NcFile *fd;

  nc::NcGroup parameter_group;
  nc::NcGroup traj_group;

  // Save dimensions
  nc::NcDim frame_dim;

  /// Variable length type for topology
  nc::NcVlenType uint_array_t;
  nc::NcVlenType double_array_t;

  /// Variable for storing time
  nc::NcVar time_var;

  /// Vlen variable for topology
  nc::NcVar topo_var;
  /// Vlen variable for coordinates
  nc::NcVar coord_var;
  /// Vlen variable for protein density
  nc::NcVar phi_var;
  /// Vlen variable for velocities
  nc::NcVar vel_var;
  /// Vlen variable for external forces
  nc::NcVar extF_var;

  /// Filepath to file
  std::string filename;
  /// Writeable status
  bool writeable;
};
} // namespace solver
} // namespace mem3dg
#endif
