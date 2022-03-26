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

#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"
#include "mem3dg/version.h"

#include <cmath>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <iostream>
#include <stdexcept>

namespace mem3dg {
namespace solver {
namespace integrator {

double Integrator::getAdaptiveCharacteristicTimeStep() {
  double currentMinimumSize = system.vpg->edgeLengths.raw().minCoeff();
  double currentMaximumForce =
      system.parameters.variation.isShapeVariation
          ? toMatrix(system.forces.mechanicalForce).cwiseAbs().maxCoeff()
          : toMatrix(system.forces.chemicalPotential).cwiseAbs().maxCoeff();

  double dt = (dt_size2_ratio * currentMinimumSize * currentMinimumSize) *
              (initialMaximumForce / currentMaximumForce);

  if (characteristicTimeStep / dt > 1e3) {
    mem3dg_runtime_message("Time step too small! May consider restarting the "
                           "simulation in small time scale");
    std::cout << "Current size / initial size = "
              << currentMinimumSize /
                     pow(characteristicTimeStep / dt_size2_ratio, 0.5)
              << std::endl;
    std::cout << "Current force / initial force = "
              << currentMaximumForce / initialMaximumForce << std::endl;
    EXIT = true;
    SUCCESS = false;
  }
  return dt;
}

double Integrator::backtrack(
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &&chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double positionProjection = 0;
  double chemicalProjection = 0;
  if (system.parameters.variation.isShapeVariation) {
    positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                          positionDirection.array())
                             .sum();
    if (positionProjection < 0) {
      std::cout << "\nBacktracking line search: positional velocity on uphill "
                   "direction, use bare "
                   "gradient! \n"
                << std::endl;
      positionDirection = toMatrix(system.forces.mechanicalForceVec);
      positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                            positionDirection.array())
                               .sum();
    }
  }
  if (system.parameters.variation.isProteinVariation) {
    chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                          chemicalDirection.array())
                             .sum();
    if (chemicalProjection < 0) {
      std::cout << "\nBacktracking line search: chemical direction on "
                   "uphill direction, "
                   "use bare "
                   "gradient! \n"
                << std::endl;
      chemicalDirection = system.forces.chemicalPotential.raw();
      chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                            chemicalDirection.array())
                               .sum();
    }
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  if (system.parameters.variation.isShapeVariation) {
    toMatrix(system.vpg->inputVertexPositions) += alpha * positionDirection;
  }
  if (system.parameters.variation.isProteinVariation) {
    system.proteinDensity.raw() += alpha * chemicalDirection;
  }
  system.time += alpha;
  system.updateConfigurations();
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (system.energy.potentialEnergy <
        (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
         c1 * alpha * (positionProjection + chemicalProjection))) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message("line search failure! Simulation "
                             "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      system.testForceComputation(alpha, toMatrix(initial_pos),
                                  toMatrix(initial_protein), previousE);
      std::cout << "\nError backtrace using characteristicTimeStep: \n"
                << std::endl;
      system.testForceComputation(characteristicTimeStep,
                                  toMatrix(system.vpg->inputVertexPositions),
                                  toMatrix(initial_protein), previousE);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    if (system.parameters.variation.isShapeVariation) {
      toMatrix(system.vpg->inputVertexPositions) =
          toMatrix(initial_pos) + alpha * positionDirection;
    }
    if (system.parameters.variation.isProteinVariation) {
      system.proteinDensity.raw() =
          toMatrix(initial_protein) + alpha * chemicalDirection;
    }
    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && ifPrintToConsole) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
    // std::cout << "mech norm: " << system.mechErrorNorm << std::endl;
    // std::cout << "chem norm: " << system.chemErrorNorm << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug)
    system.testForceComputation(alpha, toMatrix(initial_pos),
                                toMatrix(initial_protein), previousE);

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations();
  system.computePotentialEnergy();
  return alpha;
}
double Integrator::chemicalBacktrack(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &&chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double chemicalProjection = 0;
  chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                        chemicalDirection.array())
                           .sum();
  if (chemicalProjection < 0) {
    std::cout << "\nchemicalBacktracking line search: chemical direction on "
                 "uphill direction, "
                 "use bare "
                 "gradient! \n"
              << std::endl;
    chemicalDirection = system.forces.chemicalPotential.raw();
    chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                          chemicalDirection.array())
                             .sum();
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  system.proteinDensity.raw() += alpha * chemicalDirection;
  system.time += alpha;
  system.updateConfigurations();
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (system.energy.potentialEnergy <
        (previousE.potentialEnergy - c1 * alpha * chemicalProjection)) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message(
          "\nchemicalBacktrack: line search failure! Simulation "
          "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      system.testForceComputation(alpha, toMatrix(initial_pos),
                                  toMatrix(initial_protein), previousE);
      std::cout << "\nError backtrace using characteristicTimeStep: \n"
                << std::endl;
      system.testForceComputation(characteristicTimeStep, toMatrix(initial_pos),
                                  toMatrix(initial_protein), previousE);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    system.proteinDensity.raw() =
        toMatrix(initial_protein) + alpha * chemicalDirection;
    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && ifPrintToConsole) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    std::cout << "\nchemicalBacktrack: debugging \n" << std::endl;
    system.testForceComputation(alpha, toMatrix(initial_pos),
                                toMatrix(initial_protein), previousE);
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations();
  system.computePotentialEnergy();
  return alpha;
}

double Integrator::mechanicalBacktrack(
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double positionProjection = 0;
  positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                        positionDirection.array())
                           .sum();
  if (positionProjection < 0) {
    std::cout
        << "\nmechanicalBacktrack line search: positional velocity on uphill "
           "direction, use bare "
           "gradient! \n"
        << std::endl;
    positionDirection = toMatrix(system.forces.mechanicalForceVec);
    positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                          positionDirection.array())
                             .sum();
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  toMatrix(system.vpg->inputVertexPositions) += alpha * positionDirection;
  system.time += alpha;
  system.updateConfigurations();
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if ((system.energy.potentialEnergy <
         (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
          c1 * alpha * positionProjection)) &&
        std::isfinite(system.energy.potentialEnergy)) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message(
          "\nmechanicalBacktrack: line search failure! Simulation "
          "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      system.testForceComputation(alpha, toMatrix(initial_pos),
                                  toMatrix(initial_protein), previousE);
      std::cout << "\nError backtrace using characterisiticTimeStep: \n"
                << std::endl;
      system.testForceComputation(characteristicTimeStep, toMatrix(initial_pos),
                                  toMatrix(initial_protein), previousE);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    toMatrix(system.vpg->inputVertexPositions) =
        toMatrix(initial_pos) + alpha * positionDirection;

    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && ifPrintToConsole) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    std::cout << "\nmechanicalBacktrack: debugging \n" << std::endl;
    system.testForceComputation(alpha, toMatrix(initial_pos),
                                toMatrix(initial_protein), previousE);
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations();
  system.computePotentialEnergy();
  return alpha;
}

void Integrator::saveData(bool ifTrajFile, bool ifMeshFile, bool ifPrint) {
  // print in-progress information in the console
  if (ifPrint) {
    std::cout << "\n"
              << "t: " << system.time << ", "
              << "n: " << frame << ", "
              << "isSmooth: " << system.isSmooth << "\n"
              << "A, At: " << system.surfaceArea << ", "
              << system.parameters.tension.At << "/" << system.surfaceArea
              << ", "
              << "V, Vt: " << system.volume << ", "
              << system.parameters.osmotic.Vt << "/" << system.volume << ", "
              << "h: "
              << toMatrix(system.vpg->inputVertexPositions).col(2).maxCoeff()
              << "\n"
              << "E_total: " << system.energy.totalEnergy << "\n"
              << "E_kin: " << system.energy.kineticEnergy << "\n"
              << "E_pot: " << system.energy.potentialEnergy << "\n"
              << "W_ext: " << system.energy.externalWork << "\n"
              << "|e|Mech: " << system.mechErrorNorm << "\n"
              << "|e|Chem: " << system.chemErrorNorm << "\n"
              << "H: ["
              << (system.vpg->vertexMeanCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (system.vpg->vertexMeanCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]"
              << "\n"
              << "K: ["
              << (system.vpg->vertexGaussianCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (system.vpg->vertexGaussianCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]"
              << "\n"
              << "phi: [" << system.proteinDensity.raw().minCoeff() << ","
              << system.proteinDensity.raw().maxCoeff() << "]" << std::endl;

    if (EXIT)
      std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                << ", and data saved to " + outputDirectory << std::endl;
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (ifTrajFile)
    saveMutableNetcdfData();
#endif

  // save variable to richData and save ply file
  if (ifMeshFile) {
    char buffer[50];
    sprintf(buffer, ifJustGeometryPly ? "/f%d_t%d_.obj" : "/f%d_t%d_.ply",
            (int)frame, (int)system.time);
    system.saveRichData(outputDirectory + "/" + std::string(buffer),
                        ifJustGeometryPly);
  }

  frame++;
}

#ifdef MEM3DG_WITH_NETCDF
void Integrator::createMutableNetcdfFile(bool isContinue) {
  // initialize netcdf traj file
  if (isContinue) {
    mutableTrajFile.open(outputDirectory + "/" + trajFileName,
                         TrajFile::NcFile::write);
  } else {
    mutableTrajFile.createNewFile(outputDirectory + "/" + trajFileName,
                                  TrajFile::NcFile::replace);
  }
  // mutableTrajFile.writeMask(toMatrix(f.forces.forceMask).rowwise().sum());
  // if (!f.mesh->hasBoundary()) {
  //   mutableTrajFile.writeRefSurfArea(f.parameters.tension.At);
  // }
}

void Integrator::closeMutableNetcdfFile() {
  if (&mutableTrajFile != nullptr) {
    mutableTrajFile.close();
  }
}

void Integrator::saveMutableNetcdfData() {
  // scalar quantities
  // write time
  mutableTrajFile.writeTime(frame, system.time);

  // write dynamic properties
  mutableTrajFile.writeVelocity(frame, system.velocity);
  if (system.parameters.external.Kf != 0)
    mutableTrajFile.writeExternalForce(frame, system.forces.externalForceVec);

  // write static properties
  mutableTrajFile.writeCoords(frame, *system.vpg);
  mutableTrajFile.writeTopology(frame, *system.mesh);
  mutableTrajFile.writeProteinDensity(frame, system.proteinDensity);
  mutableTrajFile.sync();
}
#endif

} // namespace integrator
} // namespace solver
} // namespace mem3dg
