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

#include "mem3dg/solver/trajfile.h"

namespace mem3dg {
namespace solver {
namespace nc = ::netCDF;
using NcException = nc::exceptions::NcException;
using NcFile = nc::NcFile;

bool TrajFile::check_metadata() {
  // validate data
  std::string tmp;
  fd->getAtt(CONVENTIONS_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VALUE)
    mem3dg_runtime_error("NetCDF convention mismatch. This file does "
                         "not appear to be a valid Mem3DG trajectory.");

  fd->getAtt(CONVENTIONS_VERSION_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VERSION_VALUE)
    mem3dg_runtime_error(
        "Trajectory version mismatch. This file was generated with a "
        "different convention version.");

  return true;
}

// time & coordinate
void TrajFile::writeTime(const std::size_t idx, const double time) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  time_var.putVar({idx}, &time);
}

void TrajFile::writeIsSmooth(const std::size_t idx, const bool isSmooth) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  issmooth_var.putVar({idx}, &isSmooth);
}

void TrajFile::writeCoords(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
        &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  coord_var.putVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                   data.data());
}

void TrajFile::writeTopoFrame(const std::size_t idx,
                              const Eigen::Matrix<std::uint32_t, Eigen::Dynamic,
                                                  3, Eigen::RowMajor> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == npolygons_dim.getSize());

  topo_frame_var.putVar(
      {idx, 0, 0}, {1, npolygons_dim.getSize(), POLYGON_ORDER}, data.data());
}

double TrajFile::getTime(const std::size_t idx) const {
  assert(idx < nFrames());

  double time;
  time_var.getVar({idx}, &time);

  return time;
}

double TrajFile::getIsSmooth(const std::size_t idx) const {
  assert(idx < nFrames());

  bool isSmooth;
  issmooth_var.getVar({idx}, &isSmooth);

  return isSmooth;
}

EigenVectorX3dr TrajFile::getCoords(const std::size_t idx) const {
  assert(idx < nFrames());

  EigenVectorX3dr vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  coord_var.getVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                   vec.data());
  return vec;
}

// topology frame
Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
TrajFile::getTopoFrame(const std::size_t idx) const {
  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> vec(
      npolygons_dim.getSize(), POLYGON_ORDER);
  topo_frame_var.getVar(
      {idx, 0, 0}, {1, npolygons_dim.getSize(), POLYGON_ORDER}, vec.data());
  return vec;
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
EigenVectorX3dr TrajFile::getRefcoordinate() const {
  EigenVectorX3dr vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  refcoord.getVar({0, 0}, {nvertices_dim.getSize(), SPATIAL_DIMS}, vec.data());
  return vec;
}

// corner angles
void TrajFile::writeAngles(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == ncorners_dim.getSize());

  angle_var.putVar({idx, 0}, {1, ncorners_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getAngles(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(ncorners_dim.getSize(), 1);

  angle_var.getVar({idx, 0}, {1, ncorners_dim.getSize()}, vec.data());
  return vec;
}

// Mask
void TrajFile::writeMask(const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  assert(data.rows() == nvertices_dim.getSize());
  mask_var.putVar({0}, {nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1> TrajFile::getMask() const {
  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);
  mask_var.getVar({0}, {nvertices_dim.getSize()}, vec.data());
  return vec;
}

// velocity
void TrajFile::writeVelocity(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
        &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  vel_var.putVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                 data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS>
TrajFile::getVelocity(const std::size_t idx) const {
  assert(idx < nFrames());

  EigenVectorX3dr vec(nvertices_dim.getSize(), SPATIAL_DIMS);
  vel_var.getVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                 vec.data());
  return vec;
}

// protein density
void TrajFile::writeProteinDensity(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  phi_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getProteinDensity(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  phi_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// mean curvature
void TrajFile::writeMeanCurvature(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  meancurve_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getMeanCurvature(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  meancurve_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// Gaussian curvature
void TrajFile::writeGaussCurvature(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  gausscurve_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getGaussCurvature(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  gausscurve_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// spontaneous curvature
void TrajFile::writeSponCurvature(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  sponcurve_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getSponCurvature(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  sponcurve_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// mean - spon curvature
void TrajFile::writeH_H0_diff(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  H_H0_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getH_H0_diff(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  H_H0_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// external pressure
void TrajFile::writeExternalForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  externforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getExternalForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  externforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// chemical potential
void TrajFile::writeChemicalPotential(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  chempotential_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getChemicalPotential(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  chempotential_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// physical pressure
void TrajFile::writePhysicalForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  physforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getPhysicalForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  physforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// capillary pressure
void TrajFile::writeCapillaryForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  capforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getCapillaryForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  capforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// bending pressure
void TrajFile::writeBendingForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  bendforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getBendingForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  bendforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// inside pressure
void TrajFile::writeOsmoticForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  osmoticforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getOsmoticForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  osmoticforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// line tension pressure
void TrajFile::writeLineForce(
    const std::size_t idx,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");

  assert(data.rows() == nvertices_dim.getSize());

  lineforce_var.putVar({idx, 0}, {1, nvertices_dim.getSize()}, data.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
TrajFile::getLineForce(const std::size_t idx) const {
  assert(idx < nFrames());

  Eigen::Matrix<double, Eigen::Dynamic, 1> vec(nvertices_dim.getSize(), 1);

  lineforce_var.getVar({idx, 0}, {1, nvertices_dim.getSize()}, vec.data());
  return vec;
}

// bending energy
void TrajFile::writeBendEnergy(const std::size_t idx, const double bendEnergy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  bendener_var.putVar({idx}, &bendEnergy);
}

double TrajFile::getBendEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double bendEnergy;
  bendener_var.getVar({idx}, &bendEnergy);
  return bendEnergy;
}

// surface energy
void TrajFile::writeSurfEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  surfener_var.putVar({idx}, &Energy);
}

double TrajFile::getSurfEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  surfener_var.getVar({idx}, &Energy);
  return Energy;
}

// pressure energy
void TrajFile::writePressEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  pressener_var.putVar({idx}, &Energy);
}

double TrajFile::getPressEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  pressener_var.getVar({idx}, &Energy);
  return Energy;
}

// kinetic energy
void TrajFile::writeKineEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  kineener_var.putVar({idx}, &Energy);
}

double TrajFile::getKineEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  kineener_var.getVar({idx}, &Energy);
  return Energy;
}

// chemical energy
void TrajFile::writeAdspEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  adspener_var.putVar({idx}, &Energy);
}

double TrajFile::getAdspEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  adspener_var.getVar({idx}, &Energy);
  return Energy;
}

// line tension energy
void TrajFile::writeLineEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  lineener_var.putVar({idx}, &Energy);
}

double TrajFile::getLineEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  lineener_var.getVar({idx}, &Energy);
  return Energy;
}

// total energy
void TrajFile::writeTotalEnergy(const std::size_t idx, const double Energy) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  totalener_var.putVar({idx}, &Energy);
}

double TrajFile::getTotalEnergy(const std::size_t idx) const {
  assert(idx < nFrames());

  double Energy;
  totalener_var.getVar({idx}, &Energy);
  return Energy;
}

//  chem error norm
void TrajFile::writeChemErrorNorm(const std::size_t idx,
                                  const double ChemErrorNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  chemerrornorm_var.putVar({idx}, &ChemErrorNorm);
}

double TrajFile::getChemErrorNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double ChemErrorNorm;
  chemerrornorm_var.getVar({idx}, &ChemErrorNorm);
  return ChemErrorNorm;
}

//  error norm
void TrajFile::writeErrorNorm(const std::size_t idx, const double ErrorNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  errornorm_var.putVar({idx}, &ErrorNorm);
}

double TrajFile::getErrorNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double ErrorNorm;
  errornorm_var.getVar({idx}, &ErrorNorm);
  return ErrorNorm;
}

//  bending pressure norm
void TrajFile::writeBendNorm(const std::size_t idx, const double BendNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  bendnorm_var.putVar({idx}, &BendNorm);
}

double TrajFile::getBendNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double BendNorm;
  bendnorm_var.getVar({idx}, &BendNorm);
  return BendNorm;
}

//  capillary pressure norm
void TrajFile::writeSurfNorm(const std::size_t idx, const double ErrorNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  surfnorm_var.putVar({idx}, &ErrorNorm);
}

double TrajFile::getSurfNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double ErrorNorm;
  surfnorm_var.getVar({idx}, &ErrorNorm);
  return ErrorNorm;
}

//  inside pressure norm
void TrajFile::writePressNorm(const std::size_t idx, const double ErrorNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  pressnorm_var.putVar({idx}, &ErrorNorm);
}

double TrajFile::getPressNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double ErrorNorm;
  pressnorm_var.getVar({idx}, &ErrorNorm);
  return ErrorNorm;
}

//  line capillary pressure norm
void TrajFile::writeLineNorm(const std::size_t idx, const double ErrorNorm) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  linenorm_var.putVar({idx}, &ErrorNorm);
}

double TrajFile::getLineNorm(const std::size_t idx) const {
  assert(idx < nFrames());

  double ErrorNorm;
  linenorm_var.getVar({idx}, &ErrorNorm);
  return ErrorNorm;
}

// volume
void TrajFile::writeVolume(const std::size_t idx, const double volume) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  volume_var.putVar({idx}, &volume);
}

double TrajFile::getVolume(const std::size_t idx) const {
  assert(idx < nFrames());

  double volume;
  volume_var.getVar({idx}, &volume);
  return volume;
}

// height
void TrajFile::writeHeight(const std::size_t idx, const double height) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  height_var.putVar({idx}, &height);
}

double TrajFile::getHeight(const std::size_t idx) const {
  assert(idx < nFrames());

  double height;
  height_var.getVar({idx}, &height);
  return height;
}

// surface area
void TrajFile::writeSurfArea(const std::size_t idx, const double surfArea) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  surfarea_var.putVar({idx}, &surfArea);
}

double TrajFile::getSurfArea(const std::size_t idx) const {
  assert(idx < nFrames());

  double surfArea;
  surfarea_var.getVar({idx}, &surfArea);
  return surfArea;
}

// reference volume
void TrajFile::writeRefVolume(const double data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  refvolume.putVar({0}, &data);
}

double TrajFile::getRefVolume() const {
  double value;
  refvolume.getVar({0}, &value);
  return value;
}

// reference surface area
void TrajFile::writeRefSurfArea(const double data) {
  if (!writeable)
    mem3dg_runtime_error("Cannot write to read only file.");
  refsurfarea.putVar({0}, &data);
}

double TrajFile::getRefSurfArea() const {
  double value;
  refsurfarea.getVar({0}, &value);
  return value;
}

} // namespace solver
} // namespace mem3dg

#endif
