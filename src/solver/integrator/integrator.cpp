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
  double currentMinimumSize = system.geometry.vpg->edgeLengths.raw().minCoeff();
  double currentMaximumForce =
      system.parameters.variation.isShapeVariation
          ? system.forces.mechanicalForce.raw().cwiseAbs().maxCoeff()
          : system.forces.chemicalPotential.raw().cwiseAbs().maxCoeff();

  double dt = (dt_size2_ratio * currentMinimumSize * currentMinimumSize) *
              (initialMaximumForce / currentMaximumForce);

  if (characteristicTimeStep / dt > 1e3) {
    mem3dg_runtime_warning("Time step too small! May consider restarting the "
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

// std::tuple<EigenVectorX3dr, EigenVectorX1d> Integrator::fixPointIteration(
//     std::function<std::tuple<EigenVectorX3dr, EigenVectorX1d>(
//         EigenVectorX3dr &, EigenVectorX1d &, double)>
//         flowMap,
//     const double h, const double tolereance) {
//   EigenVectorX3dr position = toMatrix(system.geometry.vpg->vertexPositions);
//   EigenVectorX1d protein = system.proteinDensity.raw();
//   EigenVectorX3dr posRHS;
//   EigenVectorX1d proRHS;
//   std::tie<EigenVectorX3dr, EigenVectorX1d>(posRHS, proRHS) =
//       flowMap(position, protein, h);
//   double posTol = (posRHS - position).norm();
//   double proTol = (proRHS - protein).norm();
//   const double initPosTol = posTol, initProTol = proTol;
//   {
//     position = posRHS;
//     protein = proRHS;
//     std::tie<EigenVectorX3dr, EigenVectorX1d>(posRHS, proRHS) =
//         flowMap(position, protein, h);
//     posTol = (posRHS - position).norm();
//     proTol = (proRHS - protein).norm();
//   }
//   while (posTol > tolerance * initPosTol || proTol > tolerance * initProTol)
//     ;

//   return std::make_tuple(position, protein);
// }

double Integrator::backtrack(
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  // codespell: previousE
  const Energy previousE = system.energy;

  // validate the directions
  double positionProjection = 0;
  double chemicalProjection = 0;
  if (system.parameters.variation.isShapeVariation) {
    positionProjection = ((toMatrix(system.forces.conservativeForceVec) +
                           toMatrix(system.forces.externalForceVec))
                              .array() *
                          positionDirection.array())
                             .sum();
    if (positionProjection < 0)
      mem3dg_runtime_warning("Velocity on energy "
                             "uphill direction!");
  }
  if (system.parameters.variation.isProteinVariation) {
    chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                          chemicalDirection.array())
                             .sum();
    if (chemicalProjection < 0)
      mem3dg_runtime_warning("chemical evolution on energy "
                             "uphill direction!");
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.geometry.mesh);
  initial_pos = system.geometry.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.geometry.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  if (system.parameters.variation.isShapeVariation) {
    toMatrix(system.geometry.vpg->inputVertexPositions) +=
        alpha * positionDirection;
  }
  if (system.parameters.variation.isProteinVariation) {
    system.proteinDensity.raw() += alpha * chemicalDirection;
  }
  system.time += alpha;
  system.updateConfigurations();
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (system.energy.potentialEnergy <=
        (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
         c1 * alpha * (positionProjection + chemicalProjection))) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      std::cout << "\n(time=" << system.time
                << ") line search failure! Simulation "
                   "stopped."
                << std::endl;
      system.backtraceEnergyGrowth(alpha, previousE);
      // recover the initial configuration
      system.time = init_time;
      system.proteinDensity = initial_protein;
      system.geometry.vpg->inputVertexPositions = initial_pos;
      system.testConservativeForcing(alpha);
      system.testConservativeForcing(characteristicTimeStep);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    if (system.parameters.variation.isShapeVariation) {
      toMatrix(system.geometry.vpg->inputVertexPositions) =
          toMatrix(initial_pos) + alpha * positionDirection;
    }
    if (system.parameters.variation.isProteinVariation) {
      system.proteinDensity.raw() =
          initial_protein.raw() + alpha * chemicalDirection;
    }
    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.geometry.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations();
  system.computePotentialEnergy();
  return alpha;
}
double Integrator::chemicalBacktrack(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double chemicalProjection = 0;
  chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                        chemicalDirection.array())
                           .sum();
  if (chemicalProjection < 0) {
    mem3dg_runtime_warning("chemical evolution on energy "
                           "uphill direction!");
  }
  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.geometry.mesh);
  initial_pos = system.geometry.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.geometry.mesh,
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
    if (system.energy.potentialEnergy <=
        (previousE.potentialEnergy - c1 * alpha * chemicalProjection)) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      std::cout << "\n(time=" << system.time
                << ") chemicalBacktrack: line search failure! Simulation "
                   "stopped."
                << std::endl;
      system.backtraceEnergyGrowth(alpha, previousE);
      // recover the initial configuration
      system.time = init_time;
      system.proteinDensity = initial_protein;
      system.geometry.vpg->inputVertexPositions = initial_pos;
      system.testConservativeForcing(alpha);
      system.testConservativeForcing(characteristicTimeStep);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    system.proteinDensity.raw() =
        initial_protein.raw() + alpha * chemicalDirection;
    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.geometry.vpg->inputVertexPositions = initial_pos;
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
  positionProjection = ((toMatrix(system.forces.conservativeForceVec) +
                         toMatrix(system.forces.externalForceVec))
                            .array() *
                        positionDirection.array())
                           .sum();
  if (positionProjection < 0)
    mem3dg_runtime_warning("Velocity on energy "
                           "uphill direction!");

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.geometry.mesh);
  initial_pos = system.geometry.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.geometry.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  toMatrix(system.geometry.vpg->inputVertexPositions) +=
      alpha * positionDirection;
  system.time += alpha;
  system.updateConfigurations();
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if ((system.energy.potentialEnergy <=
         (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
          c1 * alpha * positionProjection)) &&
        std::isfinite(system.energy.potentialEnergy)) {
      break;
    }

    // limit of backtracking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      std::cout << "\n(time=" << system.time
                << ") mechanicalBacktrack: line search failure! Simulation "
                   "stopped."
                << std::endl;
      system.backtraceEnergyGrowth(alpha, previousE);
      // recover the initial configuration
      system.time = init_time;
      system.proteinDensity = initial_protein;
      system.geometry.vpg->inputVertexPositions = initial_pos;
      system.testConservativeForcing(alpha);
      system.testConservativeForcing(characteristicTimeStep);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    toMatrix(system.geometry.vpg->inputVertexPositions) =
        toMatrix(initial_pos) + alpha * positionDirection;

    system.time = init_time + alpha;
    system.updateConfigurations();
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.geometry.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations();
  system.computePotentialEnergy();
  return alpha;
}

void Integrator::saveData(bool ifOutputTrajFile, bool ifOutputMeshFile,
                          bool ifPrintToConsole) {
  // print in-progress information in the console
  if (ifPrintToConsole) {
    std::cout
        << "\n"
        << "t: " << system.time << ", "
        << "n: " << frame << ", "
        << "isSmooth: " << system.isSmooth << "\n"
        << "A, tension: " << system.geometry.surfaceArea << ", "
        << system.forces.surfaceTension << ", "
        << "V, pressure: " << system.geometry.volume << ", "
        << system.forces.osmoticPressure << ", "
        << "h: "
        << toMatrix(system.geometry.vpg->inputVertexPositions).col(2).maxCoeff()
        << "\n"
        << "nFaces: " << system.geometry.mesh->nFaces() << "\n"
        << "minE: " << system.geometry.vpg->edgeLengths.raw().minCoeff()
        << "; maxE: " << system.geometry.vpg->edgeLengths.raw().maxCoeff()
        << "\n"
        << "minA: " << system.geometry.vpg->faceAreas.raw().minCoeff()
        << "; maxA: " << system.geometry.vpg->faceAreas.raw().maxCoeff() << "\n"
        << "E_total: " << system.energy.totalEnergy << "\n"
        << "E_kin: " << system.energy.kineticEnergy << "\n"
        << "E_pot: " << system.energy.potentialEnergy << "\n"
        << "W_ext: " << system.energy.externalWork << "\n"
        << "|e|Mech: " << system.mechErrorNorm << "\n"
        << "|e|Chem: " << system.chemErrorNorm << "\n"
        << "H: ["
        << (system.geometry.vpg->vertexMeanCurvatures.raw().array() /
            system.geometry.vpg->vertexDualAreas.raw().array())
               .minCoeff()
        << ","
        << (system.geometry.vpg->vertexMeanCurvatures.raw().array() /
            system.geometry.vpg->vertexDualAreas.raw().array())
               .maxCoeff()
        << "]"
        << "\n"
        << "K: ["
        << (system.geometry.vpg->vertexGaussianCurvatures.raw().array() /
            system.geometry.vpg->vertexDualAreas.raw().array())
               .minCoeff()
        << ","
        << (system.geometry.vpg->vertexGaussianCurvatures.raw().array() /
            system.geometry.vpg->vertexDualAreas.raw().array())
               .maxCoeff()
        << "]"
        << "\n"
        << "phi: [" << system.proteinDensity.raw().minCoeff() << ","
        << system.proteinDensity.raw().maxCoeff() << "]"
        << "\n"
        << "sum_phi: "
        << (system.proteinDensity * system.geometry.vpg->vertexDualAreas)
               .raw()
               .sum()
        << "\n"
        << "H0: [" << system.H0.raw().minCoeff() << ","
        << system.H0.raw().maxCoeff() << "]" << std::endl;

    // report the backtracking if verbose
    if (timeStep != characteristicTimeStep && ifPrintToConsole) {
      std::cout << "timeStep: " << characteristicTimeStep << " -> " << timeStep
                << std::endl;
      system.testConservativeForcing(characteristicTimeStep);
    }

    if (EXIT)
      std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                << ", and data saved to " + outputDirectory << std::endl;
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (ifOutputTrajFile)
    saveMutableNetcdfData();
#endif

  // save variable to richData and save ply file
  if (ifOutputMeshFile) {

    std::stringstream ss;
    ss << "/" << (int)frame << "_" << (int)system.time << "_.";
    if (ifJustGeometryPly) {
      ss << "obj";
    } else {
      ss << "ply";
    }
    system.saveRichData(outputDirectory + "/" + ss.str(), ifJustGeometryPly);
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
  isMutableNetcdfFileCreated = true;
  // mutableTrajFile.writeMask(toMatrix(f.forces.forceMask).rowwise().sum());
  // if (!f.mesh->hasBoundary()) {
  //   mutableTrajFile.writeRefSurfArea(f.parameters.tension.At);
  // }
}

void Integrator::closeMutableNetcdfFile() {
  if (isMutableNetcdfFileCreated) {
    mutableTrajFile.close();
  }
}

void Integrator::saveMutableNetcdfData() {
  // scalar quantities
  // write time
  mutableTrajFile.writeTime(frame, system.time);

  // write dynamic properties
  mutableTrajFile.writeVelocity(frame, system.velocity);
  if (system.parameters.external.form != NULL)
    mutableTrajFile.writeExternalForce(frame, system.forces.externalForceVec);

  // write static properties
  mutableTrajFile.writeCoords(frame, *system.geometry.vpg);
  mutableTrajFile.writeRefCoords(frame, *system.geometry.refVpg);
  mutableTrajFile.writeTopology(frame, *system.geometry.mesh);
  mutableTrajFile.writeProteinDensity(frame, system.proteinDensity);
  mutableTrajFile.writeNotableVertex(frame, system.geometry.notableVertex);
  mutableTrajFile.sync();
}
#endif

} // namespace integrator
} // namespace solver
} // namespace mem3dg
