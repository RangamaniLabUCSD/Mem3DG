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
  for (;;) {

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData();
    }

    // Process mesh every tProcessMesh period
    if (system.time - lastProcessMesh > processMeshPeriod) {
      lastProcessMesh = system.time;
      system.mutateMesh();
      system.updateConfigurations(false);
    }

    // update geodesics every tUpdateGeodesics period
    if (system.time - lastUpdateGeodesics > updateGeodesicsPeriod) {
      lastUpdateGeodesics = system.time;
      system.updateConfigurations(true);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    if (system.time == lastProcessMesh || system.time == lastUpdateGeodesics) {
      system.time += 1e-10 * characteristicTimeStep;
      countCG = 0;
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
  if (system.parameters.external.Kf != 0) {
    mem3dg_runtime_error(
        "External force can not be applied using energy optimization")
  }
  // if (f.O.isVertexShift) {
  //   mem3dg_runtime_error(
  //       "Vertex shift is not supported for CG integration!");
  // }
}

void ConjugateGradient::status() {
  auto physicalForce = toMatrix(system.forces.mechanicalForce);

  // compute summerized forces
  system.computePhysicalForcing(timeStep);

  // compute the area contraint error
  areaDifference = abs(system.surfaceArea / system.parameters.tension.At - 1);
  if (system.parameters.osmotic.isPreferredVolume) {
    volumeDifference = abs(system.volume / system.parameters.osmotic.Vt - 1);
    reducedVolumeThreshold(EXIT, isAugmentedLagrangian, areaDifference,
                           volumeDifference, constraintTolerance, 1.3);
  } else {
    volumeDifference = abs(system.parameters.osmotic.n / system.volume /
                               system.parameters.osmotic.cam -
                           1.0);
    pressureConstraintThreshold(EXIT, isAugmentedLagrangian, areaDifference,
                                constraintTolerance, 1.3);
  }

  // exit if reached time
  if (system.time > totalTime) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
    SUCCESS = false;
  }

  // compute the free energy of the system
  system.computeTotalEnergy();

  // backtracing for error
  finitenessErrorBacktrace();
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
    system.proteinVelocity =
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
    system.proteinVelocity *= currentNormSquared / pastNormSquared;
    system.proteinVelocity +=
        system.parameters.proteinMobility * system.forces.chemicalPotential;
    pastNormSquared = currentNormSquared;
    countCG++;
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    updateAdaptiveCharacteristicStep();
  }

  // time stepping on vertex position
  if (isBacktrack) {
    timeStep = backtrack(toMatrix(system.velocity),
                         toMatrix(system.proteinVelocity), rho, c1);
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
