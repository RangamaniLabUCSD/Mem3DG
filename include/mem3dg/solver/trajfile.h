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
 * @file  trajfile.h
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

#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/typetraits.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;

// DIMENSIONS NAMES

static const std::string POLYGON_ORDER_NAME = "polygon_dims";
/// Number of vertices per polygon
static const std::size_t POLYGON_ORDER = 3;

static const std::string NPOLYGONS_NAME = "npolygons";

static const std::string SPATIAL_DIMS_NAME = "spatial";
/// Number of spatial dimensions
static const std::size_t SPATIAL_DIMS = 3;
/// Name of frames
static const std::string FRAME_NAME = "frame";
/// nvertices
static const std::string NVERTICES_NAME = "nvertices";

/// Name of conventions
static const std::string CONVENTIONS_NAME = "Conventions";
/// Convention value
static const std::string CONVENTIONS_VALUE = "Mem3DG";
/// Conventions version
static const std::string CONVENTIONS_VERSION_NAME = "ConventionsVersion";
/// Conventions version value
static const std::string CONVENTIONS_VERSION_VALUE = "0.0.1";

/// Name of the units labels
static const std::string UNITS = "units";
/// Value for length units
static const std::string LEN_UNITS = "micrometers";
/// Value for time units
static const std::string TIME_UNITS = "picoseconds";

// Data/Variable block names
/// Name of time data
static const std::string TIME_VAR = "time";
/// Name of coordinates data
static const std::string COORD_VAR = "coordinates";
/// Name of the mesh topology data
static const std::string TOPO_VAR = "topology";
/// Name of the refMesh coordinates data
static const std::string REFCOORD_VAR = "refcoordinates";
/// Name of the velocity data
static const std::string VEL_VAR = "velocities";
/// Name of the mean curvature data
static const std::string MEANCURVE_VAR = "meancurvature";
/// Name of the spontaneous curvature data
static const std::string SPONCURVE_VAR = "sponcurvature";
/// Name of the external pressure data
static const std::string EXTERNPRESS_VAR = "externpressure";
/// Name of the physical pressure data
static const std::string PHYSPRESS_VAR = "physpressure";
/// Name of the capillary pressure data
static const std::string CAPPRESS_VAR = "cappressure";
/// Name of the bending pressure data
static const std::string BENDPRESS_VAR = "bendpressure";

/**
 * @class TrajFile
 * @brief Trajectory interface to help with manipulating trajectories
 *
 */
class DLL_PUBLIC TrajFile {
private:
  TrajFile *operator=(const TrajFile &rhs) = delete;
  TrajFile(const TrajFile &rhs) = delete;

public:
  using NcFile = nc::NcFile;
  using NcException = nc::exceptions::NcException;
  using EigenVector =
      Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>;

  // template <typename T>
  // using EigenVector_T = Eigen::Matrix<T, Eigen::Dynamic, SPATIAL_DIMS,
  // Eigen::RowMajor>;

  TrajFile() = delete;

  TrajFile(TrajFile &&rhs) = default;

  /**
   * @brief Open a new file and populate it with the convention
   *
   * @param filename  Filename to save to
   * @param mesh      Mesh of interest
   * @param replace   Whether to replace an existing file or exit
   *
   * @exception netCDF::exceptions::NcExist File already exists and
   * replace/overwrite flag is not specified.
   *
   * @return TrajFile helper object to manipulate the bound NetCDF file.
   */
  static TrajFile newFile(const std::string &filename, gcs::SurfaceMesh &mesh,
                          gcs::VertexPositionGeometry &refVpg,
                          bool replace = true) {
    if (replace)
      return TrajFile(filename, mesh, refVpg, NcFile::replace);
    else
      return TrajFile(filename, mesh, refVpg, NcFile::newFile);
  };

  /**
   * @brief Open an existing NetCDF file in read/write mode
   *
   * @param filename  Filename of interest
   *
   * @exception std::runtime_error if file does not conform to the convention
   * @exception netCDF::exceptions::* If file does not exist
   *
   * @return TrajFile helper object to manipulate the bound NetCDF file.
   */
  static TrajFile openRW(const std::string &filename) {
    return TrajFile(filename, NcFile::write);
  }

  /**
   * @brief Open an existing NetCDF file in read only mode
   *
   * @param filename  Filename of interest
   *
   * @exception std::runtime_error if file does not conform to the convention
   * @exception netCDF::exceptions::* If file does not exist
   *
   * @return TrajFile helper object to manipulate the bound NetCDF file.
   */
  static TrajFile openReadOnly(const std::string &filename) {
    return TrajFile(filename, NcFile::read);
  };

  /**
   * @brief Destructor frees the bound NcFile
   */
  ~TrajFile() { delete fd; };

  /**
   * @brief Check if the bound NcFile was opened in write mode
   *
   * @return True if writable
   */
  inline bool isWriteable() { return writeable; };

  inline std::size_t getNextFrameIndex() const { return frame_dim.getSize(); };

  /**
   * @brief Validate whether or not the metadata follows convention
   *
   * @return True if okay. False otherwise.
   */
  bool check_metadata();

  void writeTime(const std::size_t idx, const double time);

  void writeCoords(const std::size_t idx, const EigenVector &data);

  std::tuple<double, EigenVector> getTimeAndCoords(const std::size_t idx) const;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
  getTopology() const;

  Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
  getRefcoordinate() const;

  void writeVelocity(const std::size_t idx, const EigenVector &data);

  Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS>
  getVelocity(const std::size_t idx) const;

  void writeMeanCurvature(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getMeanCurvature(const std::size_t idx) const;

  void writeSponCurvature(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getSponCurvature(const std::size_t idx) const;

  void
  writeExternalPressure(const std::size_t idx,
                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getExternalPressure(const std::size_t idx) const;

  void
  writePhysicalPressure(const std::size_t idx,
                        const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getPhysicalPressure(const std::size_t idx) const;

  void
  writeCapillaryPressure(const std::size_t idx,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getCapillaryPressure(const std::size_t idx) const;

  void
  writeBendingPressure(const std::size_t idx,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getBendingPressure(const std::size_t idx) const;

private:
  /**
   * @brief Private constructor for opening an existing file.
   *
   * The metadata is checked for consistency with the convention. This
   * function should only be called with fMode set to `read` or `write`.
   *
   * @param filename Path to file of interest
   * @param fMode    Mode to open file with
   */
  TrajFile(const std::string &filename, const NcFile::FileMode fMode)
      : filename(filename), // fd(new NcFile(filename, fMode)),
        writeable(fMode != NcFile::read) {

    fd = new NcFile(filename, fMode);
    check_metadata();

    frame_dim = fd->getDim(FRAME_NAME);
    npolygons_dim = fd->getDim(NPOLYGONS_NAME);
    nvertices_dim = fd->getDim(NVERTICES_NAME);
    spatial_dim = fd->getDim(SPATIAL_DIMS_NAME);
    polygon_order_dim = fd->getDim(POLYGON_ORDER_NAME);

    topology = fd->getVar(TOPO_VAR);
    refcoord = fd->getVar(REFCOORD_VAR);
    time_var = fd->getVar(TIME_VAR);
    coord_var = fd->getVar(COORD_VAR);
    vel_var = fd->getVar(VEL_VAR);
    meancurve_var = fd->getVar(MEANCURVE_VAR);
    sponcurve_var = fd->getVar(SPONCURVE_VAR);
    externpress_var = fd->getVar(EXTERNPRESS_VAR);
    physpress_var = fd->getVar(PHYSPRESS_VAR);
    cappress_var = fd->getVar(CAPPRESS_VAR);
    bendpress_var = fd->getVar(BENDPRESS_VAR);
    vel_var = fd->getVar(VEL_VAR);
  }

  /**
   * @brief Private constructor for creating a new file
   *
   * This function applies the convention. It also sets
   *
   * @param filename Path to file of interest
   * @param mesh     Mesh to store
   * @param fMode    Mode to create file with (replace, newFile)
   */
  TrajFile(const std::string &filename, gcs::SurfaceMesh &mesh,
           gcs::VertexPositionGeometry &refVpg, const NcFile::FileMode fMode)
      : filename(filename), // fd(new NcFile(filename, fMode)),
        writeable(true) {

    fd = new NcFile(filename, fMode);
    // initialize data
    fd->putAtt(CONVENTIONS_NAME, CONVENTIONS_VALUE);
    fd->putAtt(CONVENTIONS_VERSION_NAME, CONVENTIONS_VERSION_VALUE);

    frame_dim = fd->addDim(FRAME_NAME);
    npolygons_dim = fd->addDim(NPOLYGONS_NAME, mesh.nFaces());
    nvertices_dim = fd->addDim(NVERTICES_NAME, mesh.nVertices());
    spatial_dim = fd->addDim(SPATIAL_DIMS_NAME, SPATIAL_DIMS);
    polygon_order_dim = fd->addDim(POLYGON_ORDER_NAME, POLYGON_ORDER);

    // Initialize topology data block
    topology =
        fd->addVar(TOPO_VAR, nc::ncUint, {npolygons_dim, polygon_order_dim});

    // Populate topology data
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
        faceMatrix = getFaceVertexMatrix(mesh);
    std::uint32_t *topodata = faceMatrix.data();
    topology.putVar(topodata);

    // Initialize reference coordinate data block
    refcoord = fd->addVar(REFCOORD_VAR, netCDF::ncDouble,
                          {nvertices_dim, spatial_dim});

    // Populate reference coordinate data
    double *refcoorddata;
    refcoorddata = gc::EigenMap<double, 3>(refVpg.inputVertexPositions).data();
    refcoord.putVar(refcoorddata);

    time_var = fd->addVar(TIME_VAR, netCDF::ncDouble, {frame_dim});
    time_var.putAtt(UNITS, TIME_UNITS);

    coord_var = fd->addVar(COORD_VAR, netCDF::ncDouble,
                           {frame_dim, nvertices_dim, spatial_dim});
    coord_var.putAtt(UNITS, LEN_UNITS);

    vel_var = fd->addVar(VEL_VAR, netCDF::ncDouble,
                         {frame_dim, nvertices_dim, spatial_dim});

    meancurve_var =
        fd->addVar(MEANCURVE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});

    sponcurve_var =
        fd->addVar(SPONCURVE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});

    externpress_var = fd->addVar(EXTERNPRESS_VAR, netCDF::ncDouble,
                                 {frame_dim, nvertices_dim});

    physpress_var =
        fd->addVar(PHYSPRESS_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});

    cappress_var =
        fd->addVar(CAPPRESS_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});

    bendpress_var =
        fd->addVar(BENDPRESS_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
  }

  /// Bound NcFile
  NcFile *fd;

  // Save dimensions
  nc::NcDim frame_dim;
  nc::NcDim npolygons_dim;
  nc::NcDim nvertices_dim;
  nc::NcDim spatial_dim;
  nc::NcDim polygon_order_dim;

  // Save variables
  nc::NcVar topology;
  nc::NcVar refcoord;
  nc::NcVar time_var;
  nc::NcVar coord_var;
  nc::NcVar meancurve_var;
  nc::NcVar sponcurve_var;
  nc::NcVar externpress_var;
  nc::NcVar physpress_var;
  nc::NcVar cappress_var;
  nc::NcVar bendpress_var;
  nc::NcVar vel_var;

  /// Filepath to file
  std::string filename;
  /// Writeable status
  bool writeable;
};
} // namespace ddgsolver
#endif
