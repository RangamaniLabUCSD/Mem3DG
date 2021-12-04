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
#include <csignal>
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/integrator/velocity_verlet.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool VelocityVerlet::integrate() {
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
    static double lastSave;
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData();
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    march();
  }

  // return if physical simulation is sucessful
  if (!SUCCESS) {
    markFileName("_failed");
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

void VelocityVerlet::checkParameters() {
  system.meshProcessor.meshMutator.summarizeStatus();
  if (system.meshProcessor.isMeshMutate) {
    mem3dg_runtime_error(
        "Mesh mutations are currently not supported for Velocity Verlet!");
  }
}

void VelocityVerlet::status() {
  // recompute cached values
  system.updateConfigurations();

  // alias vpg quantities, which should follow the update
  auto physicalForceVec = toMatrix(system.forces.mechanicalForceVec);
  auto physicalForce = toMatrix(system.forces.mechanicalForce);
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(system.vpg->vertexNormals);

  // compute summerized forces
  getForces();

  // Compute total pressure
  // newTotalPressure = rowwiseScalarProduct((physicalForce +
  // DPDForce).array()
  // /
  //                                       f.vpg->vertexDualAreas.raw().array(),
  //                                   vertexAngleNormal_e);
  newTotalPressure =
      rowwiseScalarProduct(dpdForce.array() /
                               system.vpg->vertexDualAreas.raw().array(),
                           vertexAngleNormal_e) +
      (physicalForceVec.array().colwise() /
       system.vpg->vertexDualAreas.raw().array())
          .matrix();

  // compute the area contraint error
  areaDifference =
      (system.parameters.tension.Ksg != 0)
          ? abs(system.surfaceArea / system.parameters.tension.At - 1)
          : 0.0;

  if (system.parameters.osmotic.isPreferredVolume) {
    // compute volume constraint error
    volumeDifference =
        (system.parameters.osmotic.Kv != 0)
            ? abs(system.volume / system.parameters.osmotic.Vt - 1)
            : 0.0;
  } else {
    // compute pressure constraint error
    volumeDifference = (!system.mesh->hasBoundary())
                           ? abs(system.parameters.osmotic.n / system.volume /
                                     system.parameters.osmotic.cam -
                                 1.0)
                           : 1.0;
  }

  // exit if under error tol
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    std::cout << "\nError norm smaller than tol." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (system.time > totalTime) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
  }

  // compute the free energy of the system
  if (system.parameters.external.Kf != 0)
    system.computeExternalWork(system.time, timeStep);
  system.computeTotalEnergy();

  // backtracking for error
  finitenessErrorBacktrack();
  if (system.energy.totalEnergy > 1.05 * initialTotalEnergy) {
    std::cout
        << "\nVelocity Verlet: increasing system energy, simulation stopped!"
        << std::endl;
    SUCCESS = false;
  }
}

void VelocityVerlet::march() {
  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    updateAdaptiveCharacteristicStep();
  }

  timeStep = characteristicTimeStep;

  double hdt = 0.5 * timeStep, hdt2 = hdt * timeStep;

  // time stepping on vertex position
  toMatrix(system.vpg->inputVertexPositions) +=
      toMatrix(system.velocity) * timeStep + hdt2 * totalPressure;
  toMatrix(system.velocity) += (totalPressure + newTotalPressure) * hdt;
  totalPressure = newTotalPressure;
  system.time += timeStep;

  // time stepping on protein density
  if (system.parameters.variation.isProteinVariation) {
    system.proteinVelocity =
        system.parameters.proteinMobility * system.forces.chemicalPotential;
    system.proteinDensity += system.proteinVelocity * timeStep;
  }

  // process the mesh with regularization or mutation
  system.mutateMesh();
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
