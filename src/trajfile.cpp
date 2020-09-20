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

#ifdef MEM3DG_WITH_NETCDF

#include <cassert>
#include <iostream>
#include <netcdf>
#include <vector>

#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/trajfile.h"

namespace ddgsolver {

namespace nc = ::netCDF;
using NcException = nc::exceptions::NcException;
using NcFile = nc::NcFile;

bool TrajFile::check_metadata() {
  // validate data
  std::string tmp;
  fd->getAtt(CONVENTIONS_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VALUE)
    throw std::runtime_error("NetCDF convention mismatch. This file does "
                             "not appear to be a valid Mem3DG trajectory.");

  fd->getAtt(CONVENTIONS_VERSION_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VERSION_VALUE)
    throw std::runtime_error(
        "Trajectory version mismatch. This file was generated with a "
        "different convention version.");

  return true;
}

// time & coordinate 
void TrajFile::writeTime(const std::size_t idx, const double time) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  time_var.putVar({idx}, &time);
}

void TrajFile::writeCoords(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
        &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  coord_var.putVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                   data.data());
}

std::tuple<double, TrajFile::EigenVector>
TrajFile::getTimeAndCoords(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double time;
  time_var.getVar({idx}, &time);

  EigenVector vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  coord_var.getVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                   vec.data());
  return std::tie(time, vec);
}

// topology 
Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
TrajFile::getTopology() const {
  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> vec(
      npolygons_dim.getSize(), POLYGON_ORDER);
  topology.getVar({0, 0}, {npolygons_dim.getSize(), POLYGON_ORDER}, vec.data());
  return vec;
}

// reference coordinate
TrajFile::EigenVector
TrajFile::getRefcoordinate() const {
  EigenVector vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  refcoord.getVar({0, 0}, {nvertices_dim.getSize(), SPATIAL_DIMS}, vec.data());
  return vec;
}

// velocity 
void TrajFile::writeVelocity(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
        &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  vel_var.putVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                 data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS>
TrajFile::getVelocity(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  EigenVector vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  vel_var.getVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                 vec.data());
  return vec;
}

// mean curvature
void TrajFile::writeMeanCurvature(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  meancurve_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getMeanCurvature(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  meancurve_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// spontaneous curvature 
void TrajFile::writeSponCurvature(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  sponcurve_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getSponCurvature(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  sponcurve_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// external pressure 
void TrajFile::writeExternalPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  externpress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getExternalPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  externpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()},
                               vec.data());
  return vec;
}

// physical pressure 
void TrajFile::writePhysicalPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  physpress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()},
                              data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getPhysicalPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  physpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()},
                              vec.data());
  return vec;
}

// capillary pressure 
void TrajFile::writeCapillaryPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  cappress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()},
                               data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getCapillaryPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  cappress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()},
                               vec.data());
  return vec;
}

// bending pressure 
void TrajFile::writeBendingPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  bendpress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()},
                               data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getBendingPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  bendpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()},
                               vec.data());
  return vec;
}

} // namespace ddgsolver

#endif
