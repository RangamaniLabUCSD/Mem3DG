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

#include "system.h"

namespace mem3dg {
namespace integration {

DLL_PUBLIC void stormerVerlet(System &f, double dt, double total_time,
                              double tolerance);

/*
verbosity:
  0: no output file
  1: Traj.nc
  2: All above + console printing
  3: All above + txt + .ply
*/
DLL_PUBLIC void velocityVerlet(System &f, double dt, double init_time,
                               double total_time, double tSave,
                               double tolerance, const size_t verbosity,
                               const bool isAdaptiveStep,
                               std::string outputDir);

DLL_PUBLIC void euler(System &f, double dt, double init_time, double total_time,
                      double tSave, double tolerance, const size_t verbosity,
                      std::string outputDir, const bool isBacktrack,
                      const double rho, const double c1, const bool isAdaptiveStep);

DLL_PUBLIC void conjugateGradient(System &f, double dt, double init_time,
                                  double total_time, double tSave, double tol,
                                  double ctol, const size_t verbosity,
                                  std::string outputDir, const bool isBacktrack,
                                  const double rho, const double c1,
                                  const bool isAugmentedLagrangian,
                                  const bool isAdaptiveStep,
                                  const std::string trajFileName);

DLL_PUBLIC void
feedForwardSweep(System &f, std::vector<double> H_, std::vector<double> V_,
                 double dt, double maxTime, double tSave, double tol,
                 double ctol, std::string outputDir, const bool isBacktrack,
                 const double rho, const double c1,
                 const bool isAugmentedLagrangian, const bool isAdaptiveStep);

DLL_PUBLIC
void getParameterLog(System &f, double dt, double total_time, double tolerance,
                     double tSave, std::string inputMesh,
                     std::string outputDir);

DLL_PUBLIC void
getStatusLog(std::string nameOfFile, System &f, double dt, double time,
             std::size_t frame, double areaError, double volumeError,
             double bendingError, double faceError, double bendingEnergy,
             double surfaceEnergy, double pressureEnergy, double kineticEnergy,
             double chemicalEnergy, double lineEnergy, double totalEnergy,
             double L2ErrorNorm, bool isTuftedLaplacian, bool isProtein,
             bool isVertexShift, std::string inputMesh);

void getEnergyLog(double time, double bendingEnergy, double surfaceEnergy,
                  double pressureEnergy, double kineticEnergy,
                  double chemicalEnergy, double totalEnergy,
                  std::string outputDir);

} // namespace integration
} // namespace mem3dg
