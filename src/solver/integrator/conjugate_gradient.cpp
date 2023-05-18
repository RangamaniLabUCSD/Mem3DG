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
#include <stdexcept>

#include "Eigen/src/Core/util/Constants.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool ConjugateGradient::integrate() {
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
  for (;;) {

    // Evaluate and threshold status data
    status();

    // Save files every tSave period and print some info
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData(ifOutputTrajFile, ifOutputMeshFile, ifPrintToConsole);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    if (system.updatePrescription(lastUpdateTime, timeStep)) {
      system.time += 1e-5 * timeStep;
      countCG = 0;
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
    if (tolerance == 0) {
      markFileName(filePath, "_most", ".");
    } else {
      markFileName(filePath, "_failed", ".");
    }
  }

  return SUCCESS;
}

void ConjugateGradient::checkParameters() {
  if (system.parameters.dpd.gamma != 0) {
    mem3dg_runtime_error("DPD has to be turned off for CG integration!");
  }
  if (system.parameters.proteinMobility != 1 &&
      system.parameters.proteinMobility != 0) {
    mem3dg_runtime_error("Protein mobility constant should "
                         "be set to 1 for optimization!");
  }
  if (system.parameters.damping != 0) {
    mem3dg_runtime_error("Damping to be 0 for euler integration!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      mem3dg_runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
  if (restartPeriod < 1) {
    mem3dg_runtime_error("restartNum > 0!");
  }
  if (system.parameters.external.form != NULL) {
    mem3dg_runtime_error(
        "External force can not be applied using energy optimization")
  }
  // if (f.O.isVertexShift) {
  //   mem3dg_runtime_error(
  //       "Vertex shift is not supported for CG integration!");
  // }
}

void ConjugateGradient::status() {
  auto physicalForce = system.forces.mechanicalForce.raw();

  // compute summarized forces
  system.computeConservativeForcing();
  system.addNonconservativeForcing(timeStep);

  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    // areaDifference = abs(system.surfaceArea / system.parameters.tension.At -
    // 1); volumeDifference =
    //     (system.parameters.osmotic.isPreferredVolume)
    //         ? abs(system.volume / system.parameters.osmotic.Vt - 1)
    //         : abs(system.parameters.osmotic.n / system.volume /
    //                   system.parameters.osmotic.cam -
    //               1.0);
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
  system.computeTotalEnergy();

  // check finiteness
  if (!std::isfinite(timeStep) || !system.checkFiniteness()) {
    EXIT = true;
    SUCCESS = false;
    if (!std::isfinite(timeStep))
      mem3dg_runtime_warning("time step is not finite!");
  }
}

void ConjugateGradient::march() {
  // determine conjugate gradient direction, restart after nVertices() cycles
  if (countCG % restartPeriod == 0) {
    pastNormSquared =
        (system.parameters.variation.isShapeVariation
             ? toMatrix(system.forces.mechanicalForceVec).squaredNorm()
             : 0) +
        (system.parameters.variation.isProteinVariation
             ? system.forces.chemicalPotential.raw().squaredNorm()
             : 0);
    system.velocity = system.forces.mechanicalForceVec;
    system.proteinRateOfChange =
        system.parameters.proteinMobility * system.forces.chemicalPotential;
    countCG = 1;
  } else {
    currentNormSquared =
        (system.parameters.variation.isShapeVariation
             ? toMatrix(system.forces.mechanicalForceVec).squaredNorm()
             : 0) +
        (system.parameters.variation.isProteinVariation
             ? system.forces.chemicalPotential.raw().squaredNorm()
             : 0);
    system.velocity *= currentNormSquared / pastNormSquared;
    system.velocity += system.forces.mechanicalForceVec;
    system.proteinRateOfChange *= currentNormSquared / pastNormSquared;
    system.proteinRateOfChange +=
        system.parameters.proteinMobility * system.forces.chemicalPotential;
    pastNormSquared = currentNormSquared;
    countCG++;
  }
  system.mechErrorNorm = (toMatrix(system.velocity).array() *
                          toMatrix(system.forces.mechanicalForceVec).array())
                             .sum();
  system.chemErrorNorm = (system.proteinRateOfChange.raw().array() *
                          system.forces.chemicalPotential.raw().array())
                             .sum();

  // adjust time step if adopt adaptive time step based on mesh size
  if (ifAdaptiveStep) {
    characteristicTimeStep = getAdaptiveCharacteristicTimeStep();
  }

  // time stepping on vertex position
  if (isBacktrack) {
    timeStep = backtrack(toMatrix(system.velocity),
                         system.proteinRateOfChange.raw(), rho, c1);
  } else {
    timeStep = characteristicTimeStep;
  }
  system.geometry.vpg->inputVertexPositions += system.velocity * timeStep;
  system.proteinDensity += system.proteinRateOfChange * timeStep;
  system.time += timeStep;

  // recompute cached values
  system.updateConfigurations();
}

void ConjugateGradient::enforceAugmentedLagrangianConstraints(
    double &lambdaSG, double &lambdaV, const double dA, const double dV,
    const double tol) {
  if (ifPrintToConsole)
    std::cout << "\n["
              << "lambdaSG"
              << ", "
              << "lambdaV"
              << "] = [" << lambdaSG << ", " << lambdaV << "]";
  // update coefficient
  if (dA > tol) {
    double tension, energy;
    std::tie(tension, energy) =
        system.parameters.tension.form(system.geometry.surfaceArea);
    lambdaSG += tension;
  }
  if (dV > tol) {
    double pressure, energy;
    std::tie(pressure, energy) =
        system.parameters.osmotic.form(system.geometry.volume);
    lambdaV += pressure;
  }
  if (ifPrintToConsole)
    std::cout << " -> [" << lambdaSG << ", " << lambdaV << "]" << std::endl;
}

void ConjugateGradient::enforceIncrementalPenaltyConstraints(double &Ksg,
                                                             double &Kv,
                                                             const double dA,
                                                             const double dV,
                                                             double increment) {
  if (ifPrintToConsole)
    std::cout << "\n[Ksg, Kv] = [" << Ksg << ", " << Kv << "]";
  // update coefficient
  if (dA > constraintTolerance)
    Ksg *= increment;
  if (dV > constraintTolerance)
    Kv *= increment;
  if (ifPrintToConsole)
    std::cout << " -> [" << Ksg << ", " << Kv << "]" << std::endl;
}

} // namespace integrator
} // namespace solver
} // namespace mem3dg
