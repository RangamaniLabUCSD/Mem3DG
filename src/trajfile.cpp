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

namespace mem3dg {

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
TrajFile::EigenVector TrajFile::getRefcoordinate() const {
  EigenVector vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  refcoord.getVar({0, 0}, {nvertices_dim.getSize(), SPATIAL_DIMS}, vec.data());
  return vec;
}

// corner angles
void TrajFile::writeAngles(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == ncorners_dim.getSize());

  angle_var.putVar({idx, 0}, {1, ncorners_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getAngles(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(ncorners_dim.getSize(), 1);

  angle_var.getVar({idx, 0}, {1, ncorners_dim.getSize()}, vec.data());
  return vec;
}

// Mask
void TrajFile::writeMask(const Eigen::Matrix<int, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  assert(data.rows() == nvertices_dim.getSize());
  mask_var.putVar({0}, {nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<int, Eigen::Dynamic, 1> TrajFile::getMask() const {
  Eigen::Matrix<int, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);
  mask_var.getVar({0}, {nvertices_dim.getSize()}, vec.data());
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

// protein density
void TrajFile::writeProteinDensity(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  proteinden_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getProteinDensity(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  proteinden_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
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

// mean - spon curvature
void TrajFile::writeH_H0_diff(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  H_H0_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getH_H0_diff(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  H_H0_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
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

  externpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// physical pressure
void TrajFile::writePhysicalPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  physpress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getPhysicalPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  physpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// capillary pressure
void TrajFile::writeCapillaryPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  cappress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getCapillaryPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  cappress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// bending pressure
void TrajFile::writeBendingPressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  bendpress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getBendingPressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  bendpress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// line tension pressure
void TrajFile::writeLinePressure(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  linepress_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getLinePressure(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  linepress_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// bending energy
void TrajFile::writeBendEnergy(const std::size_t idx, const double bendEnergy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  bendener_var.putVar({idx}, &bendEnergy);
}

double TrajFile::getBendEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double bendEnergy;
  bendener_var.getVar({idx}, &bendEnergy);
  return bendEnergy;
}

// surface energy
void TrajFile::writeSurfEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  surfener_var.putVar({idx}, &Energy);
}

double TrajFile::getSurfEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  surfener_var.getVar({idx}, &Energy);
  return Energy;
}

// pressure energy
void TrajFile::writePressEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  pressener_var.putVar({idx}, &Energy);
}

double TrajFile::getPressEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  pressener_var.getVar({idx}, &Energy);
  return Energy;
}

// kinetic energy
void TrajFile::writeKineEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  kineener_var.putVar({idx}, &Energy);
}

double TrajFile::getKineEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  kineener_var.getVar({idx}, &Energy);
  return Energy;
}

// chemical energy
void TrajFile::writeChemEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  chemener_var.putVar({idx}, &Energy);
}

double TrajFile::getChemEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  chemener_var.getVar({idx}, &Energy);
  return Energy;
}

// line tension energy
void TrajFile::writeLineEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  lineener_var.putVar({idx}, &Energy);
}

double TrajFile::getLineEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  lineener_var.getVar({idx}, &Energy);
  return Energy;
}

// total energy
void TrajFile::writeTotalEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  totalener_var.putVar({idx}, &Energy);
}

double TrajFile::getTotalEnergy(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double Energy;
  totalener_var.getVar({idx}, &Energy);
  return Energy;
}

// L2 error norm
void TrajFile::writeL2ErrorNorm(const std::size_t idx,
                                const double L2ErronNorm) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  l2errornorm_var.putVar({idx}, &L2ErronNorm);
}

double TrajFile::getL2ErrorNorm(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double L2ErronNorm;
  l2errornorm_var.getVar({idx}, &L2ErronNorm);
  return L2ErronNorm;
}

// volume
void TrajFile::writeVolume(const std::size_t idx,
                                const double volume) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  l2errornorm_var.putVar({idx}, &volume);
}

double TrajFile::getVolume(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double volume;
  volume_var.getVar({idx}, &volume);
  return volume;
}

// surface area 
void TrajFile::writeSurfArea(const std::size_t idx,
                                const double surfArea) {
  if (!writeable)
    throw std::runtime_error("Cannot write to read only file.");
  surfarea_var.putVar({idx}, &surfArea);
}

double TrajFile::getSurfArea(const std::size_t idx) const {
  assert(idx < getNextFrameIndex());

  double surfArea;
  surfarea_var.getVar({idx}, &surfArea);
  return surfArea;
}

} // namespace mem3dg

#endif
