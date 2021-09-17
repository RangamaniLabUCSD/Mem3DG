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

  TrajFile() : writeable(false), fd(nullptr){};

  void open(const std::string &filename, const NcFile::FileMode fMode) {
    if ((fd != nullptr) && (fMode != NcFile::read)) {
      mem3dg_runtime_error("Cannot open an already opened ...");
    }

    fd = new NcFile(filename, fMode);
    writeable = fMode != NcFile::read;
    check_metadata();

    frame_dim = fd->getDim(FRAME_NAME);
    npolygons_dim = fd->getDim(NPOLYGONS_NAME);
    nvertices_dim = fd->getDim(NVERTICES_NAME);
    ncorners_dim = fd->getDim(NCORNERS_NAME);
    spatial_dim = fd->getDim(SPATIAL_DIMS_NAME);
    polygon_order_dim = fd->getDim(POLYGON_ORDER_NAME);

    topology = fd->getVar(TOPO_VAR);
    refcoord = fd->getVar(REFCOORD_VAR);
    angle_var = fd->getVar(ANGLE_VAR);
    time_var = fd->getVar(TIME_VAR);
    issmooth_var = fd->getVar(ISSMOOTH_VAR);
    coord_var = fd->getVar(COORD_VAR);
    topo_frame_var = fd->getVar(TOPO_FRAME_VAR);
    vel_var = fd->getVar(VEL_VAR);
    externforce_var = fd->getVar(EXTERNFORCE_VAR);
    meancurve_var = fd->getVar(MEANCURVE_VAR);
    gausscurve_var = fd->getVar(GAUSSCURVE_VAR);
    phi_var = fd->getVar(PHI_VAR);
    sponcurve_var = fd->getVar(SPONCURVE_VAR);
    chempotential_var = fd->getVar(CHEMPOTENTIAL_VAR);
    physforce_var = fd->getVar(PHYSFORCE_VAR);
    capforce_var = fd->getVar(CAPFORCE_VAR);
    bendforce_var = fd->getVar(BENDFORCE_VAR);
    osmoticforce_var = fd->getVar(OSMOTICFORCE_VAR);
    lineforce_var = fd->getVar(LINEFORCE_VAR);
    bendener_var = fd->getVar(BENDENER_VAR);
    chemerrornorm_var = fd->getVar(CHEMERRORNORM_VAR);
    errornorm_var = fd->getVar(ERRORNORM_VAR);
    bendnorm_var = fd->getVar(BENDNORM_VAR);
    surfnorm_var = fd->getVar(SURFNORM_VAR);
    pressnorm_var = fd->getVar(PRESSNORM_VAR);
    linenorm_var = fd->getVar(LINENORM_VAR);
    volume_var = fd->getVar(VOLUME_VAR);
    height_var = fd->getVar(HEIGHT_VAR);
    surfarea_var = fd->getVar(SURFAREA_VAR);
    refvolume = fd->getVar(REFVOLUME_VAR);
    refsurfarea = fd->getVar(REFSURFAREA_VAR);
    mask_var = fd->getVar(MASK_VAR);
    H_H0_var = fd->getVar(H_H0_VAR);
  }

  void createNewFile(const std::string &filename, gcs::SurfaceMesh &mesh,
                     gcs::VertexPositionGeometry &refVpg,
                     const NcFile::FileMode fMode) {
    if (fd != nullptr) {
      mem3dg_runtime_error("Cannot open an already open ...");
    }

    writeable = true;

    fd = new NcFile(filename, fMode);
    // initialize data
    fd->putAtt(CONVENTIONS_NAME, CONVENTIONS_VALUE);
    fd->putAtt(CONVENTIONS_VERSION_NAME, CONVENTIONS_VERSION_VALUE);

    frame_dim = fd->addDim(FRAME_NAME);
    npolygons_dim = fd->addDim(NPOLYGONS_NAME, mesh.nFaces());
    nvertices_dim = fd->addDim(NVERTICES_NAME, mesh.nVertices());
    ncorners_dim = fd->addDim(NCORNERS_NAME, mesh.nCorners());
    spatial_dim = fd->addDim(SPATIAL_DIMS_NAME, SPATIAL_DIMS);
    polygon_order_dim = fd->addDim(POLYGON_ORDER_NAME, POLYGON_ORDER);

    // Initialize topology data block
    topology =
        fd->addVar(TOPO_VAR, nc::ncUint, {npolygons_dim, polygon_order_dim});

    // Populate topology data
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> faceMatrix{
        mesh.getFaceVertexMatrix<std::uint32_t>()};
    std::uint32_t *topodata = faceMatrix.data();
    topology.putVar(topodata);

    // Initialize reference coordinate data block
    refcoord = fd->addVar(REFCOORD_VAR, netCDF::ncDouble,
                          {nvertices_dim, spatial_dim});
    refcoord.putAtt(UNITS, LEN_UNITS);

    // Populate reference coordinate data
    double *refcoorddata;
    refcoorddata = gc::EigenMap<double, 3>(refVpg.inputVertexPositions).data();
    refcoord.putVar(refcoorddata);

    mask_var = fd->addVar(MASK_VAR, netCDF::ncDouble, {nvertices_dim});

    time_var = fd->addVar(TIME_VAR, netCDF::ncDouble, {frame_dim});
    time_var.putAtt(UNITS, TIME_UNITS);

    issmooth_var = fd->addVar(ISSMOOTH_VAR, netCDF::ncByte, {frame_dim});

    coord_var = fd->addVar(COORD_VAR, netCDF::ncDouble,
                           {frame_dim, nvertices_dim, spatial_dim});
    coord_var.putAtt(UNITS, LEN_UNITS);

    topo_frame_var = fd->addVar(TOPO_FRAME_VAR, netCDF::ncUint,
                                {frame_dim, npolygons_dim, polygon_order_dim});

    angle_var =
        fd->addVar(ANGLE_VAR, netCDF::ncDouble, {frame_dim, ncorners_dim});

    vel_var = fd->addVar(VEL_VAR, netCDF::ncDouble,
                         {frame_dim, nvertices_dim, spatial_dim});
    vel_var.putAtt(UNITS, LEN_UNITS + TIME_UNITS + "^(-1)");

    phi_var = fd->addVar(PHI_VAR, netCDF::ncDouble,
                                {frame_dim, nvertices_dim});
    phi_var.putAtt(UNITS, LEN_UNITS + "^(-2)");

    meancurve_var =
        fd->addVar(MEANCURVE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    meancurve_var.putAtt(UNITS, LEN_UNITS + "^(-1)");

    gausscurve_var = fd->addVar(GAUSSCURVE_VAR, netCDF::ncDouble,
                                {frame_dim, nvertices_dim});
    gausscurve_var.putAtt(UNITS, LEN_UNITS + "^(-2)");

    sponcurve_var =
        fd->addVar(SPONCURVE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    sponcurve_var.putAtt(UNITS, LEN_UNITS + "^(-1)");

    externforce_var = fd->addVar(EXTERNFORCE_VAR, netCDF::ncDouble,
                                 {frame_dim, nvertices_dim});
    externforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    chempotential_var = fd->addVar(CHEMPOTENTIAL_VAR, netCDF::ncDouble,
                                   {frame_dim, nvertices_dim});
    chempotential_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    physforce_var =
        fd->addVar(PHYSFORCE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    physforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    capforce_var =
        fd->addVar(CAPFORCE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    capforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    bendforce_var =
        fd->addVar(BENDFORCE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    bendforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    osmoticforce_var = fd->addVar(OSMOTICFORCE_VAR, netCDF::ncDouble,
                                  {frame_dim, nvertices_dim});
    osmoticforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    lineforce_var =
        fd->addVar(LINEFORCE_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    lineforce_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    bendener_var = fd->addVar(BENDENER_VAR, netCDF::ncDouble, {frame_dim});
    bendener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    surfener_var = fd->addVar(SURFENER_VAR, netCDF::ncDouble, {frame_dim});
    surfener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    pressener_var = fd->addVar(PRESSENER_VAR, netCDF::ncDouble, {frame_dim});
    pressener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    kineener_var = fd->addVar(KINEENER_VAR, netCDF::ncDouble, {frame_dim});
    kineener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    adspener_var = fd->addVar(ADSPENER_VAR, netCDF::ncDouble, {frame_dim});
    adspener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    lineener_var = fd->addVar(LINEENER_VAR, netCDF::ncDouble, {frame_dim});
    lineener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    totalener_var = fd->addVar(TOTALENER_VAR, netCDF::ncDouble, {frame_dim});
    totalener_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS);

    chemerrornorm_var =
        fd->addVar(CHEMERRORNORM_VAR, netCDF::ncDouble, {frame_dim});
    chemerrornorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-1)");

    errornorm_var = fd->addVar(ERRORNORM_VAR, netCDF::ncDouble, {frame_dim});
    errornorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    bendnorm_var = fd->addVar(BENDNORM_VAR, netCDF::ncDouble, {frame_dim});
    bendnorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    surfnorm_var = fd->addVar(SURFNORM_VAR, netCDF::ncDouble, {frame_dim});
    surfnorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    pressnorm_var = fd->addVar(PRESSNORM_VAR, netCDF::ncDouble, {frame_dim});
    pressnorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    linenorm_var = fd->addVar(LINENORM_VAR, netCDF::ncDouble, {frame_dim});
    linenorm_var.putAtt(UNITS, FORCE_UNITS + LEN_UNITS + "^(-2)");

    volume_var = fd->addVar(VOLUME_VAR, netCDF::ncDouble, {frame_dim});
    volume_var.putAtt(UNITS, LEN_UNITS + "^3");

    height_var = fd->addVar(HEIGHT_VAR, netCDF::ncDouble, {frame_dim});
    height_var.putAtt(UNITS, LEN_UNITS);

    surfarea_var = fd->addVar(SURFAREA_VAR, netCDF::ncDouble, {frame_dim});
    surfarea_var.putAtt(UNITS, LEN_UNITS + "^2");

    refvolume = fd->addVar(REFVOLUME_VAR, netCDF::ncDouble);
    refvolume.putAtt(UNITS, LEN_UNITS + "^(3)");

    refsurfarea = fd->addVar(REFSURFAREA_VAR, netCDF::ncDouble);
    refsurfarea.putAtt(UNITS, LEN_UNITS + "^(2)");

    H_H0_var =
        fd->addVar(H_H0_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim});
    H_H0_var.putAtt(UNITS, LEN_UNITS + "^(-1)");
  }

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
                          bool replace = false) {
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

  inline std::size_t nFrames() const { return frame_dim.getSize(); };

  /**
   * @brief Validate whether or not the metadata follows convention
   *
   * @return True if okay. False otherwise.
   */
  bool check_metadata();

  void writeTime(const std::size_t idx, const double time);

  void writeIsSmooth(const std::size_t idx, const bool isSmooth);

  void writeCoords(const std::size_t idx, const EigenVectorX3dr &data);

  void writeTopoFrame(const std::size_t idx,
                      const Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3,
                                          Eigen::RowMajor> &data);

  double getTime(const std::size_t idx) const;

  double getIsSmooth(const std::size_t idx) const;

  EigenVectorX3dr getCoords(const std::size_t idx) const;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
  getTopoFrame(const std::size_t idx) const;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
  getTopology() const;

  Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
  getRefcoordinate() const;

  void writeAngles(const std::size_t idx,
                   const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getAngles(const std::size_t idx) const;

  void writeMask(const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1> getMask() const;

  void writeVelocity(const std::size_t idx, const EigenVectorX3dr &data);

  Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS>
  getVelocity(const std::size_t idx) const;

  void
  writeProteinDensity(const std::size_t idx,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getProteinDensity(const std::size_t idx) const;

  void writeMeanCurvature(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getMeanCurvature(const std::size_t idx) const;

  void
  writeGaussCurvature(const std::size_t idx,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getGaussCurvature(const std::size_t idx) const;

  void writeSponCurvature(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getSponCurvature(const std::size_t idx) const;

  void writeExternalForce(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getExternalForce(const std::size_t idx) const;

  void
  writeChemicalPotential(const std::size_t idx,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getChemicalPotential(const std::size_t idx) const;

  void writePhysicalForce(const std::size_t idx,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getPhysicalForce(const std::size_t idx) const;

  void
  writeCapillaryForce(const std::size_t idx,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getCapillaryForce(const std::size_t idx) const;

  void writeBendingForce(const std::size_t idx,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getBendingForce(const std::size_t idx) const;

  void writeOsmoticForce(const std::size_t idx,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getOsmoticForce(const std::size_t idx) const;

  void writeLineForce(const std::size_t idx,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getLineForce(const std::size_t idx) const;

  void writeH_H0_diff(const std::size_t idx,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &data);

  Eigen::Matrix<double, Eigen::Dynamic, 1>
  getH_H0_diff(const std::size_t idx) const;

  void writeBendEnergy(const std::size_t idx, const double bendEnergy);

  double getBendEnergy(const std::size_t idx) const;

  void writeSurfEnergy(const std::size_t idx, const double Energy);

  double getSurfEnergy(const std::size_t idx) const;

  void writePressEnergy(const std::size_t idx, const double Energy);

  double getPressEnergy(const std::size_t idx) const;

  void writeKineEnergy(const std::size_t idx, const double Energy);

  double getKineEnergy(const std::size_t idx) const;

  void writeAdspEnergy(const std::size_t idx, const double Energy);

  double getAdspEnergy(const std::size_t idx) const;

  void writeLineEnergy(const std::size_t idx, const double Energy);

  double getLineEnergy(const std::size_t idx) const;

  void writeTotalEnergy(const std::size_t idx, const double Energy);

  double getTotalEnergy(const std::size_t idx) const;

  void writeChemErrorNorm(const std::size_t idx, const double ChemErrorNorm);

  double getChemErrorNorm(const std::size_t idx) const;

  void writeErrorNorm(const std::size_t idx, const double ErrorNorm);

  double getErrorNorm(const std::size_t idx) const;

  void writeBendNorm(const std::size_t idx, const double BendNorm);

  double getBendNorm(const std::size_t idx) const;

  void writeSurfNorm(const std::size_t idx, const double SurfNorm);

  double getSurfNorm(const std::size_t idx) const;

  void writePressNorm(const std::size_t idx, const double PressNorm);

  double getPressNorm(const std::size_t idx) const;

  void writeLineNorm(const std::size_t idx, const double LineNorm);

  double getLineNorm(const std::size_t idx) const;

  void writeVolume(const std::size_t idx, const double volume);

  double getVolume(const std::size_t idx) const;

  void writeHeight(const std::size_t idx, const double volume);

  double getHeight(const std::size_t idx) const;

  void writeSurfArea(const std::size_t idx, const double surfArea);

  double getSurfArea(const std::size_t idx) const;

  void writeRefVolume(const double data);

  double getRefVolume() const;

  void writeRefSurfArea(const double data);

  double getRefSurfArea() const;

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
    open(filename, fMode);
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
    createNewFile(filename, mesh, refVpg, fMode);
  }

  /// Bound NcFile
  NcFile *fd;

  // Save dimensions
  nc::NcDim frame_dim;
  nc::NcDim npolygons_dim;
  nc::NcDim nvertices_dim;
  nc::NcDim ncorners_dim;
  nc::NcDim spatial_dim;
  nc::NcDim polygon_order_dim;

  // Save variables
  nc::NcVar topology;
  nc::NcVar refcoord;
  nc::NcVar time_var;
  nc::NcVar issmooth_var;
  nc::NcVar coord_var;
  nc::NcVar topo_frame_var;
  nc::NcVar angle_var;
  nc::NcVar vel_var;
  nc::NcVar phi_var;
  nc::NcVar meancurve_var;
  nc::NcVar gausscurve_var;
  nc::NcVar sponcurve_var;
  nc::NcVar externforce_var;
  nc::NcVar chempotential_var;
  nc::NcVar physforce_var;
  nc::NcVar capforce_var;
  nc::NcVar bendforce_var;
  nc::NcVar osmoticforce_var;
  nc::NcVar lineforce_var;
  nc::NcVar bendener_var;
  nc::NcVar surfener_var;
  nc::NcVar pressener_var;
  nc::NcVar kineener_var;
  nc::NcVar adspener_var;
  nc::NcVar lineener_var;
  nc::NcVar totalener_var;
  nc::NcVar chemerrornorm_var;
  nc::NcVar errornorm_var;
  nc::NcVar bendnorm_var;
  nc::NcVar surfnorm_var;
  nc::NcVar pressnorm_var;
  nc::NcVar linenorm_var;
  nc::NcVar volume_var;
  nc::NcVar height_var;
  nc::NcVar surfarea_var;
  nc::NcVar refvolume;
  nc::NcVar refsurfarea;
  nc::NcVar mask_var;
  nc::NcVar H_H0_var;

  /// Filepath to file
  std::string filename;
  /// Writeable status
  bool writeable;
};
} // namespace solver
} // namespace mem3dg
#endif
