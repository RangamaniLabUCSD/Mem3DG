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

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/bfgs.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool BFGS::integrate() {

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
    std::cout << "Initialized integrator and the output trajactory is "
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

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    march();
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

void BFGS::checkParameters() {
  if (system.parameters.dpd.gamma != 0) {
    mem3dg_runtime_error("DPD has to be turned off for BFGS integration!");
  }
  // if (f.O.isVertexShift) {
  //   mem3dg_runtime_error(
  //       "Vertex shift is not supported for BFGS integration!");
  // }
  if (!isBacktrack) {
    mem3dg_runtime_error("Backtracking is required for BFGS integration");
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
  if (system.parameters.external.Kf != 0) {
    mem3dg_runtime_error(
        "External force can not be applied using energy optimization")
  }
}

void BFGS::status() {
  auto physicalForceVec = toMatrix(system.forces.mechanicalForceVec);
  auto physicalForce = toMatrix(system.forces.mechanicalForce);

  // compute summerized forces
  system.computePhysicalForcing(timeStep);

  // update
  if (system.time != initialTime || ifRestart) {
    EigenVectorX1d y = -flatten(physicalForceVec) + pastPhysicalForce;
    double sTy = (s.transpose() * y);
    hess_inv +=
        (s * s.transpose()) * (sTy + y.transpose() * hess_inv * y) / sTy / sTy -
        (hess_inv * y * s.transpose() + s * y.transpose() * hess_inv) / sTy;

    EigenVectorX1d y_protein =
        -system.forces.chemicalPotential.raw() + pastPhysicalForce_protein;
    double sTy_protein = (s_protein.transpose() * y_protein);
    hess_inv_protein += (s_protein * s_protein.transpose()) *
                            (sTy_protein + y_protein.transpose() *
                                               hess_inv_protein * y_protein) /
                            sTy_protein / sTy_protein -
                        (hess_inv_protein * y_protein * s_protein.transpose() +
                         s_protein * y_protein.transpose() * hess_inv_protein) /
                            sTy_protein;
  }
  pastPhysicalForce = flatten(physicalForceVec);
  pastPhysicalForce_protein = system.forces.chemicalPotential.raw();
  // std::cout << "if equal: "
  //           << (unflatten<3>(flatten(physicalForceVec)).array() ==
  //               physicalForceVec.array())
  //           << std::endl;

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
    // thresholding, exit if fulfilled and iterate if not
    reducedVolumeThreshold(EXIT, isAugmentedLagrangian, areaDifference,
                           volumeDifference, ctol, 1.3);
  } else {
    // compute pressure constraint error
    volumeDifference = (!system.mesh->hasBoundary())
                           ? abs(system.parameters.osmotic.n / system.volume /
                                     system.parameters.osmotic.cam -
                                 1.0)
                           : 1.0;
    // thresholding, exit if fulfilled and iterate if not
    pressureConstraintThreshold(EXIT, isAugmentedLagrangian, areaDifference,
                                ctol, 1.3);
  }

  // exit if reached time
  if (system.time > totalTime) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
    SUCCESS = false;
  }

  // compute the free energy of the system
  system.computeTotalEnergy();

  // backtracking for error
  finitenessErrorBacktrace();
}

void BFGS::march() {
  if (system.time == lastSave && system.time != initialTime) {
    // process the mesh with regularization or mutation
    system.mutateMesh();
    system.updateConfigurations(true);

    system.time += 1e-10 * characteristicTimeStep;
    ifRestart = true;
    hess_inv.setIdentity();
    hess_inv_protein.setIdentity();
  } else {
    // map the raw eigen datatype for computation
    auto f_velocity_e = toMatrix(system.velocity);
    auto f_forces_mechanicalForceVec_e =
        toMatrix(system.forces.mechanicalForceVec);

    f_velocity_e = unflatten<3>(
        (hess_inv * flatten(f_forces_mechanicalForceVec_e)).eval());
    toMatrix(system.proteinVelocity) =
        hess_inv_protein * system.forces.chemicalPotential.raw();

    // adjust time step if adopt adaptive time step based on mesh size and force
    // magnitude
    if (isAdaptiveStep) {
      updateAdaptiveCharacteristicStep();
    }

    // time stepping on vertex position
    timeStep =
        backtrack(f_velocity_e, toMatrix(system.proteinVelocity), rho, c1);
    system.vpg->inputVertexPositions += system.velocity * timeStep;
    system.proteinDensity += system.proteinVelocity * timeStep;
    system.time += timeStep;
    s = timeStep * flatten(f_velocity_e);
    s_protein = timeStep * toMatrix(system.proteinVelocity);

    // regularization
    if (system.meshProcessor.isMeshRegularize) {
      system.computeRegularizationForce();
      system.vpg->inputVertexPositions.raw() +=
          system.forces.regularizationForce.raw();
    }

    // recompute cached values
    system.updateConfigurations(false);
  }
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
