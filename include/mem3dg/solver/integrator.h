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

#pragma once

#include <geometrycentral/surface/geometry.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>

#include "force.h"

namespace ddgsolver {
namespace integration {

DLL_PUBLIC std::tuple<double, double, double, double, double> getFreeEnergy(Force &f);

DLL_PUBLIC double
getL2ErrorNorm(Eigen::SparseMatrix<double> M,
               Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure);

DLL_PUBLIC void stormerVerlet(Force &f, double dt, double total_time,
                              double tolerance);

DLL_PUBLIC void velocityVerlet(Force &f, double dt, double total_time,
                               double tolerance, double closeZone,
                               double increment, double maxKv, double maxKsg,
                               double tSave, double tMollify,
                               std::string inputMesh, std::string outputDir);

DLL_PUBLIC void getParameterLog(Force &f, double dt,
                                double total_time,
                                double tolerance, double tSave,
                                std::string inputMesh, std::string outputDir);

DLL_PUBLIC void getStatusLog(std::string nameOfFile, Force &f, double dt, double time,
                              std::size_t frame, double areaError, double volumeError,
                              double bendingError, double faceError,
                              double bendingEnergy, double surfaceEnergy, double pressureEnergy, double chemicalEnergy,
                              double totalEnergy,
                              double L2ErrorNorm, bool isTuftedLaplacian, bool isProtein, bool isVertexShift,
                              std::string inputMesh);

} // namespace integration
} // namespace ddgsolver
