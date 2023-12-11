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

    // step forward
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
    markFileName(filePath, "_failed", ".");
  }

  return SUCCESS;
}

void VelocityVerlet::checkParameters() {
  // system.meshProcessor.meshMutator.summarizeStatus();
  // if (system.meshProcessor.isMeshMutate) {
  //   mem3dg_runtime_error(
  //       "Mesh mutations are currently not supported for Velocity Verlet!");
  // }
  if (system.parameters.damping == 0) {
    mem3dg_runtime_error("Expect nonzero damping force for Velocity Verlet "
                         "integration! Note that 0 < damping/time step < 1!");
  }
  if (isBacktrack && system.parameters.dpd.gamma != 0) {
    mem3dg_runtime_warning(
        "Fluctuation can lead to failure in backtracking algorithm!");
  }
}

void VelocityVerlet::status() {
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

  // check energy increase
  if (isCapEnergy) {
    if (system.energy.totalEnergy - system.energy.proteinInteriorPenalty >
        1.05 * initialTotalEnergy) {
      if (ifPrintToConsole)
        std::cout << "\nVelocity Verlet: increasing system energy, simulation "
                     "stopped! E_total="
                  << system.energy.totalEnergy -
                         system.energy.proteinInteriorPenalty
                  << ", E_init=" << initialTotalEnergy << " (w/o inPE)"
                  << std::endl;
      EXIT = true;
      SUCCESS = false;
    }
  }
}

void VelocityVerlet::march() {
  // compute protein velocity, which is independent of time
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

  // adjust characteristic time step if adopt adaptive time step based on mesh
  // size
  if (ifAdaptiveStep) {
    characteristicTimeStep = getAdaptiveCharacteristicTimeStep();
  }

  // backtracking to obtain stable time step or assumed characteristic time
  // step, ignore higher order acceleration effect
  if (isBacktrack && system.time > (1e-5 * characteristicTimeStep)) {
    double timeStep_mech = std::numeric_limits<double>::max(),
           timeStep_chem = std::numeric_limits<double>::max();
    if (system.parameters.variation.isShapeVariation)
      timeStep_mech = mechanicalBacktrack(toMatrix(system.velocity), rho, c1);
    if (system.parameters.variation.isProteinVariation)
      timeStep_chem =
          chemicalBacktrack(system.proteinRateOfChange.raw(), rho, c1);
    timeStep = (timeStep_chem < timeStep_mech) ? timeStep_chem : timeStep_mech;
  } else {
    timeStep = characteristicTimeStep;
  }

  // march the system with increment time step obtained above
  double hdt = 0.5 * timeStep, hdt2 = hdt * timeStep;

  // stepping on vertex position
  system.geometry.vpg->inputVertexPositions +=
      system.velocity * timeStep +
      hdt2 *
          pastMechanicalForceVec; // x_{i+1} = x_i + dt_i * v_i + 0.5 * (dt_i)^2
                                  // * a_i

  // stepping on protein density
  system.proteinDensity += system.proteinRateOfChange * timeStep;

  // stepping on time
  system.time += timeStep; // t_{i+1} = t_i + dt_i

  // velocity predictor for force calculation
  system.velocity +=
      timeStep * pastMechanicalForceVec; // v_{i+1} = v_i + a_i * dt_i

  // compute summarized forces
  system.computeConservativeForcing();
  system.addNonconservativeForcing(
      timeStep); // a_{i+1} at (x_{i+1}, v_{i+1}, dt_i)

  // stepping on velocity from velocity predictor
  system.velocity +=
      (system.forces.mechanicalForceVec - pastMechanicalForceVec) *
      hdt; // v_{i+1} = v_{i+1} + 0.5
           // * (a_{i+1} - a_i) * dt

  system.mechErrorNorm = (toMatrix(system.velocity).array() *
                          toMatrix(system.forces.mechanicalForceVec).array())
                             .sum();

  // cache current force
  pastMechanicalForceVec = system.forces.mechanicalForceVec; // a_i <-- a_{i+1}

  // recompute cached values
  system.updateConfigurations();
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
