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

  signal(SIGINT, signalHandler);

#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (verbosity > 0) {
    // createNetcdfFile();
    createMutableNetcdfFile();
    // print to console
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
        if (verbosity > 2) {
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

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info; save data before exit
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData();
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // Process mesh every tProcessMesh period
    if (system.time - lastProcessMesh > (processMeshPeriod * timeStep)) {
      lastProcessMesh = system.time;
      system.mutateMesh();
      if (system.meshProcessor.meshRegularizer.isSmoothenMesh)
        system.smoothenMesh(timeStep);
      system.updateConfigurations(false);
    }

    // update geodesics every tUpdateGeodesics period
    if (system.time - lastUpdateGeodesics >
        (updateGeodesicsPeriod * timeStep)) {
      lastUpdateGeodesics = system.time;
      system.updateConfigurations(true);
    }

    // step forward
    if (system.time == lastProcessMesh || system.time == lastUpdateGeodesics) {
      system.time += 1e-10 * characteristicTimeStep;
    } else {
      march();
    }
  }

  // return if optimization is sucessful
  if (!SUCCESS) {
    if (tolerance == 0) {
      markFileName("_most");
    } else {
      markFileName("_failed");
    }
  }
  // stop the timer and report time spent
#ifdef __linux__
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal integration time: " << duration << " seconds"
              << std::endl;
  }
#endif

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
  // compute summerized forces
  system.computePhysicalForcing(timeStep);

  // compute the contraint error
  areaDifference = abs(system.surfaceArea / system.parameters.tension.At - 1);
  volumeDifference = (system.parameters.osmotic.isPreferredVolume)
                         ? abs(system.volume / system.parameters.osmotic.Vt - 1)
                         : abs(system.parameters.osmotic.n / system.volume /
                                   system.parameters.osmotic.cam -
                               1.0);

  // exit if under error tolerance
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    std::cout << "\nError norm smaller than tolerance." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (system.time > totalTime) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
    SUCCESS = false;
  }

  // compute the free energy of the system
  if (system.parameters.external.Kf != 0)
    system.computeExternalWork(system.time, timeStep);
  system.computeTotalEnergy();

  // backtracking for error
  finitenessErrorBacktrace();
}

void Euler::march() {
  // compute force, which is equivalent to velocity
  system.velocity = system.forces.mechanicalForceVec;
  system.proteinVelocity =
      system.parameters.proteinMobility * system.forces.chemicalPotential;

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    characteristicTimeStep = updateAdaptiveCharacteristicStep();
  }

  // time stepping on vertex position
  if (isBacktrack) {
    double timeStep_mech,
        timeStep_chem = std::numeric_limits<double>::infinity();
    if (system.parameters.variation.isShapeVariation)
      timeStep_mech = mechanicalBacktrack(toMatrix(system.velocity), rho, c1);
    if (system.parameters.variation.isProteinVariation)
      timeStep_chem =
          chemicalBacktrack(toMatrix(system.proteinVelocity), rho, c1);
    timeStep = (timeStep_chem < timeStep_mech) ? timeStep_chem : timeStep_mech;
  } else {
    timeStep = characteristicTimeStep;
  }
  system.vpg->inputVertexPositions += system.velocity * timeStep;
  system.proteinDensity += system.proteinVelocity * timeStep;
  system.time += timeStep;

  // regularization
  if (system.meshProcessor.isMeshRegularize) {
    system.computeRegularizationForce();
    system.vpg->inputVertexPositions.raw() +=
        system.forces.regularizationForce.raw();
  }

  // recompute cached values
  system.updateConfigurations(false);
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
