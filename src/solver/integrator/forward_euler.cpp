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

#include <Eigen/Core>
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool Euler::integrate() {

  if (ifDisableIntegrate)
    mem3dg_runtime_error("integrate() is disabled for current construction!");

  signal(SIGINT, signalHandler);

  double initialTime = system.time, lastComputeAvoidingForce = system.time,
         lastSave = system.time;
  std::map<std::string, double> lastUpdateTime{{"geodesics", system.time},
                                               {"mutateMesh", system.time},
                                               {"protein", system.time},
                                               {"notableVertex", system.time},
                                               {"mask", system.time}};

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (ifOutputTrajFile) {
    createMutableNetcdfFile(isContinuation);
    if (ifPrintToConsole)
      std::cout << "Initialized NetCDF file at "
                << outputDirectory + "/" + trajFileName << std::endl;
  }
#endif

  // time integration loop
  const double avoidStrength = system.parameters.selfAvoidance.mu;
  for (;;) {

    // turn on/off self-avoidance; outside status-march-cycle; before savedata
    // to write selfAvoidance
    if (avoidStrength != 0) {
      if ((system.time - lastComputeAvoidingForce) >
              system.parameters.selfAvoidance.p * system.projectedCollideTime ||
          system.time - lastSave >= savePeriod || system.time == initialTime ||
          EXIT) {
        lastComputeAvoidingForce = system.time;
        system.parameters.selfAvoidance.mu = avoidStrength;
        if (ifPrintToConsole) {
          std::cout << "computing avoiding force at "
                    << "t = " << system.time << std::endl;
          std::cout << "projected collision is " << system.projectedCollideTime
                    << std::endl;
          std::cout << "time step is " << timeStep << std::endl;
        }
      } else {
        system.parameters.selfAvoidance.mu = 0;
      }
    }

    // Evaluate and threshold status data
    status();

    // Save files every tSave period and print some info; save data before exit
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData(ifOutputTrajFile, ifOutputMeshFile, ifPrintToConsole);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    if (system.updatePrescription(lastUpdateTime, timeStep)) {
      system.time += 1e-5 * timeStep;
    } else {
      march();
    }
  }

#ifdef MEM3DG_WITH_NETCDF
  if (ifOutputTrajFile) {
    closeMutableNetcdfFile();
    if (ifPrintToConsole)
      std::cout << "Closed NetCDF file" << std::endl;
  }
#endif

  // return if optimization is successful
  if (!SUCCESS && ifOutputTrajFile) {
    std::string filePath = outputDirectory;
    filePath.append("/");
    filePath.append(trajFileName);
  }

  return SUCCESS;
}

void Euler::checkParameters() {
  if (system.parameters.dpd.gamma != 0) {
    mem3dg_runtime_error("DPD has to be turned off for euler integration!");
  }
  if (system.parameters.damping != 0) {
    mem3dg_runtime_error("Damping to be 0 for euler integration!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      mem3dg_runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
}

void Euler::status() {
  // compute summearized forces
  system.computeConservativeForcing();
  system.addNonconservativeForcing(timeStep);

  // exit if under error tolerance
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    if (ifPrintToConsole)
      std::cout << "\nError norm smaller than tolerance." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (system.time > totalTime) {
    if (ifPrintToConsole)
      std::cout << "\nReached time." << std::endl;
    EXIT = true;
  }

  // compute the free energy of the system
  if (system.parameters.external.form != NULL)
    system.computeExternalWork(system.time, timeStep);
  system.computeTotalEnergy();

  // check finiteness
  if (!std::isfinite(timeStep) || !system.checkFiniteness()) {
    EXIT = true;
    SUCCESS = false;
    if (!std::isfinite(timeStep))
      mem3dg_runtime_warning("time step is not finite!");
  }
}

// std::tuple<EigenVectorX3dr, EigenVectorX1d>
// Euler::flowMap(EigenVectorX3dr &position, EigenVectorX1d &protein, double h)
// {
//   toMatrix(system.geometry.vpg->inputVertexPositions) = position;
//   system.proteinDensity.raw() = protein;

//   system.updateConfigurations();
//   system.computeConservativeForcing();
//   system.addNonconservativeForcing(h);

//   system.velocity = system.forces.mechanicalForceVec;

//   if (system.parameters.variation.isProteinVariation) {
//     if (system.parameters.variation.isProteinConservation) {
//       system.proteinRateOfChange.raw() =
//           system.parameters.proteinMobility *
//           system.geometry.vpg->hodge0Inverse *
//           system.geometry.vpg->d0.transpose() *
//           system.computeInPlaneFluxForm(system.forces.chemicalPotential.raw());
//     } else {
//       system.proteinRateOfChange = system.parameters.proteinMobility *
//                                    system.forces.chemicalPotential /
//                                    system.geometry.vpg->vertexDualAreas;
//     }
//   }

//   return std::make_tuple(
//       toMatrix(system.geometry.vpg->inputVertexPositions) +
//           toMatrix(system.velocity) * h,
//       (system.proteinDensity + system.proteinRateOfChange * h).raw());
// }

void Euler::march() {
  // compute velocity, which are independent of time

  // fixPointIteration(flowMap, characteristicTimeStep, 0.01);

  system.velocity = system.forces.mechanicalForceVec;
  system.mechErrorNorm = (toMatrix(system.velocity).array() *
                          toMatrix(system.forces.mechanicalForceVec).array())
                             .sum();

  if (system.parameters.variation.isProteinVariation) {
    if (system.parameters.variation.isProteinConservation) {
      system.proteinRateOfChange.raw() =
          system.parameters.proteinMobility *
          system.geometry.vpg->hodge0Inverse *
          system.geometry.vpg->d0.transpose() *
          system.computeInPlaneFluxForm(system.forces.chemicalPotential.raw());
    } else {
      system.proteinRateOfChange = system.parameters.proteinMobility *
                                   system.forces.chemicalPotential /
                                   system.geometry.vpg->vertexDualAreas;
    }
    system.chemErrorNorm = (system.proteinRateOfChange.raw().array() *
                            system.forces.chemicalPotential.raw().array())
                               .sum();
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (ifAdaptiveStep) {
    characteristicTimeStep = getAdaptiveCharacteristicTimeStep();
  }

  // backtracking to obtain stable time step
  if (isBacktrack) {
    double timeStep_mech = std::numeric_limits<double>::max(),
           timeStep_chem = std::numeric_limits<double>::max();
    if (system.parameters.variation.isShapeVariation)
      timeStep_mech = mechanicalBacktrack(toMatrix(system.velocity), rho, c1);
    if (system.parameters.variation.isProteinVariation)
      timeStep_chem =
          chemicalBacktrack(system.proteinRateOfChange.raw(), rho, c1);
    timeStep = std::min(timeStep_chem, timeStep_mech);
    // (timeStep_chem < timeStep_mech) ? timeStep_chem : timeStep_mech;
  } else {
    timeStep = characteristicTimeStep;
  }
  system.geometry.vpg->inputVertexPositions += system.velocity * timeStep;
  system.proteinDensity += system.proteinRateOfChange * timeStep;
  system.time += timeStep;

  // recompute cached values
  system.updateConfigurations();
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
